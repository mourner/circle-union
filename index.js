import Flatbush from 'flatbush';
import {within} from 'geoflatbush';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km

/**
 * @typedef {Object} State
 * @property {number} n
 * @property {Float64Array} lng
 * @property {Float64Array} lat
 * @property {Float64Array} r       radius in km
 * @property {Float64Array} cx
 * @property {Float64Array} cy
 * @property {Float64Array} cz      center unit vectors
 * @property {Float64Array} cosR
 * @property {Float64Array} sinR    cos/sin of angular radius ρ = r/R
 * @property {Flatbush} index
 */

/**
 * §0 setup. Reorders circles by radius (largest first) so the owner rule
 * collapses to `i < j` and engulfers are always processed before what they
 * engulf. Precomputes 3-D unit vectors + angular-radius trig (all the
 * transcendentals, once) and builds the spatial index. All downstream hot
 * loops are pure dot/cross products.
 * @param {Float64Array} lng
 * @param {Float64Array} lat
 * @param {Float64Array} r radius in km
 * @returns {State}
 */
export function build(lng, lat, r) {
    const n = lng.length;

    // sort indices by radius, descending
    const order = new Uint32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
    order.sort((a, b) => r[b] - r[a]);

    const slng = new Float64Array(n);
    const slat = new Float64Array(n);
    const sr = new Float64Array(n);
    const cx = new Float64Array(n);
    const cy = new Float64Array(n);
    const cz = new Float64Array(n);
    const cosR = new Float64Array(n);
    const sinR = new Float64Array(n);

    const index = new Flatbush(n);

    for (let i = 0; i < n; i++) {
        const o = order[i];
        const lngi = lng[o], lati = lat[o], ri = r[o];
        slng[i] = lngi;
        slat[i] = lati;
        sr[i] = ri;

        const latR = lati * RAD, lngR = lngi * RAD;
        const cosLat = Math.cos(latR);
        cx[i] = cosLat * Math.cos(lngR);
        cy[i] = cosLat * Math.sin(lngR);
        cz[i] = Math.sin(latR);

        const rho = ri / R;
        cosR[i] = Math.cos(rho);
        sinR[i] = Math.sin(rho);

        index.add(lngi, lati);
    }
    index.finish();

    return {n, lng: slng, lat: slat, r: sr, cx, cy, cz, cosR, sinR, index};
}

/**
 * §1–§3. Single radius-descending pass over all owner pairs (`i < j`):
 *   §2  classify each overlapping pair; flag any circle fully engulfed by
 *       another as `covered` (globally redundant — excluded from all further
 *       work, both its own searches and as a target of others').
 *   §3  for each *properly* intersecting pair, solve the two planes ∩ sphere
 *       for the two boundary points `p± = α·cᵢ + β·cⱼ ± γ·n` and store them once
 *       with a stable integer ID shared by both circles' arc lists.
 *
 * Point storage is interleaved `[x,y,z]` (points are always consumed as whole
 * vectors). Pair storage is `[i, j, baseId]` per proper pair; the pair's two
 * points are `baseId` (the `+γ` root) and `baseId + 1` (the `−γ` root).
 *
 * @param {State} state
 * @returns {{covered: Uint8Array, pairCount: number, coveredCount: number,
 *   points: Float64Array, pointCount: number, pairs: Int32Array}}
 */
export function scan(state) {
    const {n, lng, lat, r, cx, cy, cz, cosR, sinR, index} = state;
    const covered = new Uint8Array(n);
    let pairCount = 0, coveredCount = 0;

    // growable interleaved point buffer (x,y,z) and pair buffer (i,j,baseId)
    let points = new Float64Array(1 << 16); // 3 slots/point
    let np = 0;                             // point count
    let pairs = new Int32Array(1 << 15);    // 3 slots/pair
    let npairs = 0;

    for (let i = 0; i < n; i++) {
        if (covered[i]) continue;
        const ri = r[i];
        const xi = cx[i], yi = cy[i], zi = cz[i], cosRi = cosR[i], sinRi = sinR[i];

        within(index, lng[i], lat[i], 2 * ri, (j) => {
            if (j <= i || covered[j]) return false; // self, the larger owner, or already-dropped
            const xj = cx[j], yj = cy[j], zj = cz[j], cosRj = cosR[j], sinRj = sinR[j];
            const cij = xi * xj + yi * yj + zi * zj; // = cos(angular distance)
            const cosProd = cosRi * cosRj;
            const sinProd = sinRi * sinRj;
            if (cij <= cosProd - sinProd) return false; // disjoint: cij ≤ cos(ρi+ρj)

            if (cij >= cosProd + sinProd) {             // engulf: cij ≥ cos(ρi−ρj), i ⊇ j
                covered[j] = 1;
                coveredCount++;
                return false;
            }

            // proper intersection — solve for the two points
            const det = 1 - cij * cij;                  // = |n|² = |cᵢ×cⱼ|²
            const alpha = (cosRi - cij * cosRj) / det;
            const beta = (cosRj - cij * cosRi) / det;
            let g2 = (1 - alpha * alpha - beta * beta - 2 * alpha * beta * cij) / det;
            if (g2 < 0) g2 = 0;                          // tangency / roundoff guard
            const gamma = Math.sqrt(g2);

            // n = cᵢ × cⱼ
            const nx = yi * zj - zi * yj;
            const ny = zi * xj - xi * zj;
            const nz = xi * yj - yi * xj;

            // mid = α·cᵢ + β·cⱼ ; p± = mid ± γ·n
            const mx = alpha * xi + beta * xj;
            const my = alpha * yi + beta * yj;
            const mz = alpha * zi + beta * zj;
            const gx = gamma * nx, gy = gamma * ny, gz = gamma * nz;

            if (np * 3 + 6 > points.length) {
                const grown = new Float64Array(points.length * 2);
                grown.set(points); points = grown;
            }
            if (npairs * 3 + 3 > pairs.length) {
                const grown = new Int32Array(pairs.length * 2);
                grown.set(pairs); pairs = grown;
            }

            const baseId = np;
            let p = np * 3;
            points[p++] = mx + gx; points[p++] = my + gy; points[p++] = mz + gz; // p+ = baseId
            points[p++] = mx - gx; points[p++] = my - gy; points[p++] = mz - gz; // p- = baseId+1
            np += 2;

            let q = npairs * 3;
            pairs[q++] = i; pairs[q++] = j; pairs[q++] = baseId;
            npairs++;
            pairCount++;

            return false; // we account inline; keep within()'s result array empty
        });
    }

    return {covered, pairCount, coveredCount, points, pointCount: np, pairs};
}
