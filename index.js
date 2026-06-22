import Flatbush from 'flatbush';
import {within} from 'geoflatbush';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km
const TWO_PI = 2 * Math.PI;

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
 * @property {Float64Array} ux
 * @property {Float64Array} uy
 * @property {Float64Array} uz      local east unit vector (tangent frame)
 * @property {Float64Array} vx
 * @property {Float64Array} vy
 * @property {Float64Array} vz      local north unit vector (tangent frame)
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
    const ux = new Float64Array(n), uy = new Float64Array(n), uz = new Float64Array(n);
    const vx = new Float64Array(n), vy = new Float64Array(n), vz = new Float64Array(n);

    const index = new Flatbush(n);

    for (let i = 0; i < n; i++) {
        const o = order[i];
        const lngi = lng[o], lati = lat[o], ri = r[o];
        slng[i] = lngi;
        slat[i] = lati;
        sr[i] = ri;

        const latR = lati * RAD, lngR = lngi * RAD;
        const cosLat = Math.cos(latR), sinLat = Math.sin(latR);
        const cosLng = Math.cos(lngR), sinLng = Math.sin(lngR);
        cx[i] = cosLat * cosLng;
        cy[i] = cosLat * sinLng;
        cz[i] = sinLat;

        // local tangent frame at the center: east u, north v (both ⊥ c and each other,
        // u × v = c so increasing θ is CCW seen from outside the sphere). θ=0 → east.
        ux[i] = -sinLng;          uy[i] = cosLng;           uz[i] = 0;
        vx[i] = -sinLat * cosLng; vy[i] = -sinLat * sinLng; vz[i] = cosLat;

        const rho = ri / R;
        cosR[i] = Math.cos(rho);
        sinR[i] = Math.sin(rho);

        index.add(lngi, lati);
    }
    index.finish();

    return {n, lng: slng, lat: slat, r: sr, cx, cy, cz, cosR, sinR, ux, uy, uz, vx, vy, vz, index};
}

/**
 * §1–§3. Single radius-descending pass over all owner pairs (`i < j`):
 *   §2  classify each overlapping pair; flag any circle fully engulfed by
 *       another as `covered` (globally redundant — excluded from all further
 *       work, both its own searches and as a target of others').
 *   §3  for each *properly* intersecting pair, solve the two planes ∩ sphere
 *       for the two boundary points `p± = α·cᵢ + β·cⱼ ± γ·n` and store them once
 *       with a stable integer ID shared by both circles' arc lists, and `union`
 *       the pair in a disjoint-set forest — the overlap graph's connected
 *       components *are* the union's connected components, which §7 uses to nest
 *       holes into shells by grouping rather than by point-in-ring search.
 *
 * Point storage is interleaved `[x,y,z]` (points are always consumed as whole
 * vectors). Pair storage is `[i, j, baseId]` per proper pair; the pair's two
 * points are `baseId` (the `+γ` root) and `baseId + 1` (the `−γ` root).
 *
 * `component` holds a compacted component id (0…componentCount−1) per active
 * circle; engulfed (`covered`) circles get −1 (dropped, never contribute a ring).
 *
 * @param {State} state
 * @returns {{covered: Uint8Array, pairCount: number, coveredCount: number,
 *   points: Float64Array, pointCount: number, pairs: Int32Array,
 *   component: Int32Array, componentCount: number}}
 */
export function scan(state) {
    const {n, lng, lat, r, cx, cy, cz, cosR, sinR, index} = state;
    const covered = new Uint8Array(n);
    let pairCount = 0, coveredCount = 0;

    // disjoint-set forest over circles (union by size + path halving). Each proper
    // pair unions its two endpoints; roots partition active circles into components.
    const parent = new Int32Array(n);
    const setSize = new Int32Array(n);
    for (let i = 0; i < n; i++) { parent[i] = i; setSize[i] = 1; }
    const find = (x) => {
        while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    const union = (a, b) => {
        let ra = find(a), rb = find(b);
        if (ra === rb) return;
        if (setSize[ra] < setSize[rb]) { const t = ra; ra = rb; rb = t; }
        parent[rb] = ra; setSize[ra] += setSize[rb];
    };

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

            union(i, j); // same connected component of the union

            return false; // we account inline; keep within()'s result array empty
        });
    }

    // compact roots → dense component ids over active circles
    const component = new Int32Array(n).fill(-1);
    let componentCount = 0;
    for (let i = 0; i < n; i++) {
        if (covered[i]) continue;
        const root = find(i);
        if (component[root] === -1) component[root] = componentCount++;
        component[i] = component[root];
    }

    return {covered, pairCount, coveredCount, points, pointCount: np, pairs, component, componentCount};
}

/**
 * §4. Per-circle interval complement → boundary arcs.
 *
 * Each proper pair `(i, j)` from §3 covers a single angular interval on circle `i`
 * (the part of `∂i` inside disk `j`) and, symmetrically, one on circle `j`. We bound
 * each interval by the two shared intersection points' bearings in the circle's local
 * frame — `θ(p) = atan2(p·v, p·u)` — and pick the arc containing the bearing toward the
 * other center (deepest inside the other disk). Unioning all covered intervals on a
 * circle and taking the complement gives its surviving boundary arcs; each arc endpoint
 * is a §3 intersection point, referenced by its stable ID so §5 stitches by ID.
 *
 * The union/complement is a circular depth sweep: +1 at each covered interval's CCW
 * start, −1 at its end; boundary arcs are the maximal runs where depth returns to 0.
 * `baseDepth` seeds the sweep with the coverage straddling the atan2 seam (−π).
 *
 *   - active circle, no proper neighbors → one full-circle arc (startId = −1)
 *   - active circle fully covered by ≥2 neighbors jointly → zero arcs
 *   - engulfed (`covered`) circle → skipped (handled in §2)
 *
 * @param {State} state
 * @param {{covered: Uint8Array, pairCount: number, points: Float64Array, pairs: Int32Array}} scanResult
 */
export function arcs(state, scanResult) {
    const {n, cx, cy, cz, ux, uy, uz, vx, vy, vz} = state;
    const {covered, pairCount, points, pairs} = scanResult;

    // Covered-interval events grouped per circle: theta, delta (+1 start / −1 end), and the
    // intersection-point ID at that endpoint. Each proper pair yields one interval per circle
    // = two events; baseDepth[c] counts intervals wrapping the seam (CCW start > end).
    const E = pairCount * 4;
    const evTheta = new Float64Array(E);
    const evDelta = new Int8Array(E);
    const evId = new Int32Array(E);
    const evCircle = new Int32Array(E);
    const baseDepth = new Int32Array(n);
    let ne = 0;

    // record circle c's covered interval bounded by points pa (id idA) and pb (id idB),
    // selecting the arc toward (dirx,diry,dirz) — the bearing to the other circle's center.
    const addInterval = (c, ax, ay, az, idA, bx, by, bz, idB, dirx, diry, dirz) => {
        const eux = ux[c], euy = uy[c], euz = uz[c];
        const nvx = vx[c], nvy = vy[c], nvz = vz[c];
        const thA = Math.atan2(ax * nvx + ay * nvy + az * nvz, ax * eux + ay * euy + az * euz);
        const thB = Math.atan2(bx * nvx + by * nvy + bz * nvz, bx * eux + by * euy + bz * euz);
        const mid = Math.atan2(dirx * nvx + diry * nvy + dirz * nvz, dirx * eux + diry * euy + dirz * euz);

        // CCW from thA to thB; covered arc is the side containing `mid`
        let dab = thB - thA; if (dab < 0) dab += TWO_PI;
        let dam = mid - thA; if (dam < 0) dam += TWO_PI;
        let s, e, sId, eId;
        if (dam <= dab) { s = thA; e = thB; sId = idA; eId = idB; }
        else { s = thB; e = thA; sId = idB; eId = idA; }

        if (s > e) baseDepth[c]++; // interval wraps the seam
        evTheta[ne] = s; evDelta[ne] = 1; evId[ne] = sId; evCircle[ne] = c; ne++;
        evTheta[ne] = e; evDelta[ne] = -1; evId[ne] = eId; evCircle[ne] = c; ne++;
    };

    for (let pi = 0; pi < pairCount; pi++) {
        const q = pi * 3;
        const i = pairs[q], j = pairs[q + 1], baseId = pairs[q + 2];
        const pp = baseId * 3;
        const ax = points[pp], ay = points[pp + 1], az = points[pp + 2];     // p+ (baseId)
        const bx = points[pp + 3], by = points[pp + 4], bz = points[pp + 5]; // p− (baseId+1)
        addInterval(i, ax, ay, az, baseId, bx, by, bz, baseId + 1, cx[j], cy[j], cz[j]);
        addInterval(j, ax, ay, az, baseId, bx, by, bz, baseId + 1, cx[i], cy[i], cz[i]);
    }

    // counting sort event indices by circle → contiguous per-circle ranges
    const off = new Int32Array(n + 1);
    for (let k = 0; k < ne; k++) off[evCircle[k] + 1]++;
    for (let c = 0; c < n; c++) off[c + 1] += off[c];
    const order = new Int32Array(ne);
    const cursor = off.slice(0, n);
    for (let k = 0; k < ne; k++) order[cursor[evCircle[k]]++] = k;

    // boundary arcs: circle index, CCW [thetaStart, thetaEnd], endpoint point IDs.
    // startId = −1 marks a full-circle arc (no endpoints).
    let cap = 1 << 12;
    let arcCircle = new Int32Array(cap);
    let arcThetaStart = new Float64Array(cap);
    let arcThetaEnd = new Float64Array(cap);
    let arcStartId = new Int32Array(cap);
    let arcEndId = new Int32Array(cap);
    let arcCount = 0, fullCount = 0;

    const pushArc = (c, ts, te, sId, eId) => {
        if (arcCount === cap) {
            cap *= 2;
            arcCircle = grow32(arcCircle, cap); arcStartId = grow32(arcStartId, cap); arcEndId = grow32(arcEndId, cap);
            arcThetaStart = grow64(arcThetaStart, cap); arcThetaEnd = grow64(arcThetaEnd, cap);
        }
        arcCircle[arcCount] = c; arcThetaStart[arcCount] = ts; arcThetaEnd[arcCount] = te;
        arcStartId[arcCount] = sId; arcEndId[arcCount] = eId; arcCount++;
    };

    for (let c = 0; c < n; c++) {
        if (covered[c]) continue;
        const lo = off[c], hi = off[c + 1];

        if (lo === hi) { // active, no proper neighbors → whole circle is one arc
            pushArc(c, 0, TWO_PI, -1, -1);
            fullCount++;
            continue;
        }

        // insertion-sort this circle's event indices by theta (ties: +1 before −1)
        for (let a = lo + 1; a < hi; a++) {
            const key = order[a];
            const kt = evTheta[key], kd = evDelta[key];
            let b = a - 1;
            while (b >= lo && (evTheta[order[b]] > kt || (evTheta[order[b]] === kt && evDelta[order[b]] < kd))) {
                order[b + 1] = order[b]; b--;
            }
            order[b + 1] = key;
        }

        // sweep: boundary arcs are runs where coverage depth is 0
        let depth = baseDepth[c];
        let gapTheta = 0, gapId = -1, haveGap = false;   // a boundary arc opened (coverage hit 0)
        let seamTheta = 0, seamId = -1, haveSeam = false; // the start closing the seam-wrapping arc
        for (let a = lo; a < hi; a++) {
            const k = order[a];
            if (evDelta[k] === 1) {
                if (depth === 0) {
                    if (haveGap) { pushArc(c, gapTheta, evTheta[k], gapId, evId[k]); haveGap = false; }
                    else { seamTheta = evTheta[k]; seamId = evId[k]; haveSeam = true; }
                }
                depth++;
            } else {
                depth--;
                if (depth === 0) { gapTheta = evTheta[k]; gapId = evId[k]; haveGap = true; }
            }
        }
        if (haveGap && haveSeam) pushArc(c, gapTheta, seamTheta, gapId, seamId); // arc straddling the seam
    }

    return {arcCount, fullCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId, arcEndId};
}

/** @param {Int32Array} a @param {number} cap */
function grow32(a, cap) { const g = new Int32Array(cap); g.set(a); return g; }
/** @param {Float64Array} a @param {number} cap */
function grow64(a, cap) { const g = new Float64Array(cap); g.set(a); return g; }
