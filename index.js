import Flatbush from 'flatbush';
import {within} from 'geoflatbush';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km
const TWO_PI = 2 * Math.PI;

/**
 * @typedef {[number, number, number, number, number]} Arc
 *   `[lng, lat, radius, startAngle, endAngle]` — lng/lat°, radius km, angles in radians measured from east
 *   CCW, pre-unwrapped so `endAngle ≥ startAngle` and `sweep = endAngle − startAngle ∈ (0, 2π]`. A full
 *   circle is one `[…, 0, 2π]`.
 * @typedef {Arc[]} Ring          ordered CCW (interior on the left)
 * @typedef {Ring[]} Polygon      shell first, then holes
 * @typedef {Polygon[]} Topology  what `arcs()` returns
 */

/**
 * Union of geographic disks. A Flatbush-style builder: reserve a circle count, `add` each circle, then read
 * the result as exact arc topology (`arcs()`) or sampled GeoJSON (`geojson()`). The heavy pipeline runs once
 * on the first read and is cached; only `geojson`'s sampling step re-runs per call.
 */
export class CircleUnion {
    /** @param {number} numItems number of circles to reserve space for */
    constructor(numItems) {
        if (!(numItems >= 0)) throw new Error('numItems must be a non-negative number.');
        this._lng = new Float64Array(numItems);
        this._lat = new Float64Array(numItems);
        this._r = new Float64Array(numItems);
        this._pos = 0;
        /** @type {Topology | null} */
        this._topology = null; // cached arc topology, invalidated by `add`
    }

    /**
     * Add a circle. Returns its index. Throws past the reserved count.
     * @param {number} lng longitude in degrees
     * @param {number} lat latitude in degrees
     * @param {number} r radius in km
     * @returns {number}
     */
    add(lng, lat, r) {
        if (this._pos >= this._lng.length) throw new Error('Added more circles than reserved.');
        const i = this._pos++;
        this._lng[i] = lng; this._lat[i] = lat; this._r[i] = r;
        this._topology = null;
        return i;
    }

    /** Run the heavy pipeline once and cache the arc topology. */
    _compute() {
        if (this._topology) return;
        const n = this._pos;
        if (n === 0) { this._topology = []; return; }
        const state = build(this._lng.subarray(0, n), this._lat.subarray(0, n), this._r.subarray(0, n));
        const scanResult = scan(state);
        const arcResult = arcs(state, scanResult);
        this._topology = stitch(state, scanResult, arcResult);
    }

    /**
     * Exact arc topology — resolution-independent. `[polygon, ...]`, polygon = `[ring, ...]`, ring =
     * `[arc, ...]`, arc = `[lng, lat, radius, startAngle, endAngle]` (see the `Arc` typedef). Cached.
     * @returns {Topology}
     */
    arcs() {
        this._compute();
        return /** @type {Topology} */ (this._topology);
    }

    /**
     * GeoJSON `MultiPolygon` sampled from the arc topology.
     * @param {{tolerance?: number, minPoints?: number}} [options] `tolerance`: max arc↔chord deviation in km
     *   (default 0.005 ≈ 5 m); `minPoints`: floor on vertices per full circle (default 24).
     * @returns {{type: 'MultiPolygon', coordinates: number[][][][]}}
     */
    geojson(options) {
        this._compute();
        return sample(/** @type {Topology} */ (this._topology), options);
    }
}

/**
 * Setup. Sorts circles by radius descending so the "larger circle owns the pair" rule collapses to `i < j`
 * and any engulfer is processed before what it engulfs. Precomputes per-circle 3-D unit vectors and
 * angular-radius trig (all the transcendentals, once) and builds the spatial index, leaving every downstream
 * hot loop to pure dot/cross products.
 * @param {Float64Array} lng
 * @param {Float64Array} lat
 * @param {Float64Array} r radius in km
 */
function build(lng, lat, r) {
    const n = lng.length;

    // sort indices by radius, descending
    const order = new Uint32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
    order.sort((a, b) => r[b] - r[a]);

    const slng = new Float64Array(n);
    const slat = new Float64Array(n);
    const sr = new Float64Array(n);
    const cx = new Float64Array(n); // center unit vectors
    const cy = new Float64Array(n);
    const cz = new Float64Array(n);
    const cosR = new Float64Array(n); // cos/sin of angular radius ρ = r/R
    const sinR = new Float64Array(n);
    const ux = new Float64Array(n), uy = new Float64Array(n); // local east unit vector
    const vx = new Float64Array(n), vy = new Float64Array(n), vz = new Float64Array(n); // local north unit vector

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

        // local tangent frame at the center: east u, north v (both ⊥ c and each other, u × v = c so
        // increasing θ is CCW seen from outside the sphere). θ=0 → east. East is horizontal, so its
        // z-component is identically 0 and is never stored.
        ux[i] = -sinLng;          uy[i] = cosLng;
        vx[i] = -sinLat * cosLng; vy[i] = -sinLat * sinLng; vz[i] = cosLat;

        const rho = ri / R;
        cosR[i] = Math.cos(rho);
        sinR[i] = Math.sin(rho);

        index.add(lngi, lati);
    }
    index.finish();

    return {n, lng: slng, lat: slat, r: sr, cx, cy, cz, cosR, sinR, ux, uy, vx, vy, vz, index};
}

/**
 * Single radius-descending pass over all owner pairs (`i < j`):
 *   - classify each overlapping pair; flag any circle fully engulfed by another as `covered` (globally
 *     redundant — excluded from all further work, both its own searches and as a target of others').
 *   - for each *properly* intersecting pair, solve the two planes ∩ sphere for the two boundary points
 *     `p± = α·cᵢ + β·cⱼ ± γ·n`, store them once under a stable integer ID shared by both circles' arc lists,
 *     and `union` the pair in a disjoint-set forest — the overlap graph's connected components *are* the
 *     union's, which `stitch` uses to nest holes into shells by grouping rather than by point-in-ring search.
 *
 * Point storage is interleaved `[x,y,z]` (points are always consumed as whole vectors). Pair storage is
 * `[i, j, baseId]` per proper pair; the pair's two points are `baseId` (the `+γ` root) and `baseId + 1`
 * (the `−γ` root). `component` holds a compacted component id (0…componentCount−1) per active circle;
 * engulfed (`covered`) circles get −1 (dropped, never contribute a ring).
 *
 * @param {ReturnType<typeof build>} state
 * @returns {{covered: Uint8Array, pairCount: number, coveredCount: number,
 *   points: Float64Array, pointCount: number, pairs: Int32Array,
 *   component: Int32Array, componentCount: number}}
 */
function scan(state) {
    const {n, lng, lat, r, cx, cy, cz, cosR, sinR, index} = state;
    const covered = new Uint8Array(n);
    let pairCount = 0, coveredCount = 0;

    // disjoint-set forest over circles (union by size + path halving, see dsuFind/dsuUnion).
    // Each proper pair unions its two endpoints; roots partition active circles into components.
    const parent = new Int32Array(n);
    const setSize = new Int32Array(n);
    for (let i = 0; i < n; i++) { parent[i] = i; setSize[i] = 1; }

    // growable interleaved buffers: points (x, y, z) and pairs (i, j, baseId).
    // pairCount doubles as the pair-write cursor.
    let points = new Float64Array(1 << 16);
    let pairs = new Int32Array(1 << 15);
    let np = 0; // point count

    // Drop exact-duplicate circles up front (co-located tower rows with bit-identical lng/lat are common in
    // real data). Must precede the pair sweep: otherwise one larger circle crossing several duplicates mints
    // bit-identical intersection points under distinct IDs, desyncing the ring handoff in `stitch`. An O(n)
    // open-addressing hash set keyed on the exact centre does it — the first (largest-radius) circle at a
    // location wins, later duplicates are marked covered. Only *bit-identical* centres need this; near
    // coincidences mint distinct points and are caught by the engulf test below, so the match is exact float
    // equality (round() only derives hash bits).
    let cap = 1; while (cap < n * 2) cap <<= 1;
    const mask = cap - 1;
    const slot = new Int32Array(cap).fill(-1);
    for (let i = 0; i < n; i++) {
        let h = (Math.imul(Math.round(lat[i] * 1e7), 0x9e3779b1) ^ Math.round(lng[i] * 1e7)) & mask;
        while (slot[h] !== -1 && (lat[slot[h]] !== lat[i] || lng[slot[h]] !== lng[i])) h = (h + 1) & mask;
        if (slot[h] === -1) slot[h] = i;                    // new location → representative
        else { covered[i] = 1; coveredCount++; }            // same centre as an earlier, larger circle
    }

    // Radius-descending owner sweep. `within` returns each owner's candidate neighbors (no filterFn → it
    // just collects them), and we classify each pair inline. Neighbors come back in the index's traversal
    // order, so processing is deterministic.
    for (let i = 0; i < n; i++) {
        if (covered[i]) continue;
        const xi = cx[i], yi = cy[i], zi = cz[i], cosRi = cosR[i], sinRi = sinR[i];
        const neighbors = within(index, lng[i], lat[i], 2 * r[i]);

        for (let t = 0; t < neighbors.length; t++) {
            const j = neighbors[t];
            if (j <= i || covered[j]) continue; // self, the larger owner, or already-dropped
            const xj = cx[j], yj = cy[j], zj = cz[j], cosRj = cosR[j], sinRj = sinR[j];
            const cij = xi * xj + yi * yj + zi * zj; // = cos(angular distance)
            const cosProd = cosRi * cosRj;
            const sinProd = sinRi * sinRj;
            if (cij <= cosProd - sinProd) continue; // disjoint: cij ≤ cos(ρi+ρj)

            if (cij >= cosProd + sinProd) {         // engulf: cij ≥ cos(ρi−ρj), i ⊇ j
                covered[j] = 1;
                coveredCount++;
                continue;
            }

            // proper intersection — solve p± = α·cᵢ + β·cⱼ ± γ·n for the two points, n = cᵢ × cⱼ. |n|² = sin²d
            // comes from the cross product directly, NOT from 1−cij²: for close or small circles cij≈1 and 1−cij²
            // cancels catastrophically (loses ~8 digits, throwing points hundreds of m off and desyncing the
            // bearings `arcs` derives from them). The α/β numerators are likewise written with 1−cosd = 2sin²(d/2)
            // to dodge the same cancellation in cosRᵢ−cij·cosRⱼ. Keeps points on-circle to ~1e-15.
            const nx = yi * zj - zi * yj;
            const ny = zi * xj - xi * zj;
            const nz = xi * yj - yi * xj;
            const sin2d = nx * nx + ny * ny + nz * nz;   // = |cᵢ×cⱼ|² = sin²(angular dist)
            const dAng = Math.atan2(Math.sqrt(sin2d), cij);
            const sh = Math.sin(dAng / 2);
            const oneMinusCosd = 2 * sh * sh;
            const alpha = ((cosRi - cosRj) + cosRj * oneMinusCosd) / sin2d;
            const beta = ((cosRj - cosRi) + cosRi * oneMinusCosd) / sin2d;
            let g2 = (1 - alpha * alpha - beta * beta - 2 * alpha * beta * cij) / sin2d;
            if (g2 < 0) g2 = 0;                          // tangency / roundoff guard
            const gamma = Math.sqrt(g2);

            // mid = α·cᵢ + β·cⱼ ; p± = mid ± γ·n
            const mx = alpha * xi + beta * xj;
            const my = alpha * yi + beta * yj;
            const mz = alpha * zi + beta * zj;
            const gx = gamma * nx, gy = gamma * ny, gz = gamma * nz;

            if (np * 3 + 6 > points.length) points = grow64(points, points.length * 2);
            if (pairCount * 3 + 3 > pairs.length) pairs = grow32(pairs, pairs.length * 2);

            // store the pair's two boundary points (p+ = baseId, p− = baseId + 1) and the pair
            const baseId = np, p = np * 3;
            points[p] = mx + gx; points[p + 1] = my + gy; points[p + 2] = mz + gz;
            points[p + 3] = mx - gx; points[p + 4] = my - gy; points[p + 5] = mz - gz;
            np += 2;

            const q = pairCount * 3;
            pairs[q] = i; pairs[q + 1] = j; pairs[q + 2] = baseId;
            pairCount++;

            dsuUnion(parent, setSize, i, j); // i and j fall in the same connected component
        }
    }

    // compact roots → dense component ids over active circles
    const component = new Int32Array(n).fill(-1);
    let componentCount = 0;
    for (let i = 0; i < n; i++) {
        if (covered[i]) continue;
        const root = dsuFind(parent, i);
        if (component[root] === -1) component[root] = componentCount++;
        component[i] = component[root];
    }

    return {covered, pairCount, coveredCount, points, pointCount: np, pairs, component, componentCount};
}

/**
 * Per-circle interval complement → boundary arcs.
 *
 * Each proper pair `(i, j)` from `scan` covers a single angular interval on circle `i` (the part of `∂i`
 * inside disk `j`) and, symmetrically, one on circle `j`. We bound each interval by the two shared
 * intersection points' bearings in the circle's local frame — `θ(p) = atan2(p·v, p·u)` — and pick the arc
 * containing the bearing toward the other center (deepest inside the other disk). Unioning all covered
 * intervals on a circle and taking the complement gives its surviving boundary arcs; each arc endpoint is an
 * intersection point, referenced by its stable ID so `stitch` joins arcs by ID.
 *
 * The union/complement is a circular depth sweep: +1 at each covered interval's CCW start, −1 at its end;
 * boundary arcs are the maximal runs where depth returns to 0. `baseDepth` seeds the sweep with the coverage
 * straddling the atan2 seam (−π).
 *
 *   - active circle, no proper neighbors → one full-circle arc (startId = −1)
 *   - active circle fully covered by ≥2 neighbors jointly → zero arcs
 *   - engulfed (`covered`) circle → skipped
 *
 * @param {ReturnType<typeof build>} state
 * @param {{covered: Uint8Array, pairCount: number, points: Float64Array, pairs: Int32Array}} scanResult
 */
function arcs(state, scanResult) {
    const {n, cx, cy, cz, ux, uy, vx, vy, vz} = state;
    const {covered, pairCount, points, pairs} = scanResult;

    // Covered-interval events grouped per circle: theta, delta (+1 start / −1 end), and the intersection-point
    // ID at that endpoint. Each proper pair yields one interval per circle = two events; baseDepth[c] counts
    // intervals wrapping the seam (CCW start > end).
    const E = pairCount * 4;
    const evTheta = new Float64Array(E);
    const evDelta = new Int8Array(E);
    const evId = new Int32Array(E);
    const evCircle = new Int32Array(E);
    const baseDepth = new Int32Array(n);
    let ne = 0;

    // record circle c's covered interval bounded by points pa (id idA) and pb (id idB), selecting the arc
    // toward (dirx,diry,dirz) — the bearing to the other circle's center.
    /** @param {number} c
     *  @param {number} ax @param {number} ay @param {number} az @param {number} idA
     *  @param {number} bx @param {number} by @param {number} bz @param {number} idB
     *  @param {number} dirx @param {number} diry @param {number} dirz */
    const addInterval = (c, ax, ay, az, idA, bx, by, bz, idB, dirx, diry, dirz) => {
        const eux = ux[c], euy = uy[c];                  // east, z-component ≡ 0
        const nvx = vx[c], nvy = vy[c], nvz = vz[c];
        const thA = Math.atan2(ax * nvx + ay * nvy + az * nvz, ax * eux + ay * euy);
        const thB = Math.atan2(bx * nvx + by * nvy + bz * nvz, bx * eux + by * euy);
        const mid = Math.atan2(dirx * nvx + diry * nvy + dirz * nvz, dirx * eux + diry * euy);

        // CCW from thA to thB; covered arc is the side containing `mid`
        let dab = thB - thA; if (dab < 0) dab += TWO_PI;
        let dam = mid - thA; if (dam < 0) dam += TWO_PI;
        let s, e, sId, eId;
        if (dam <= dab) {
            s = thA; e = thB; sId = idA; eId = idB;
        } else {
            s = thB; e = thA; sId = idB; eId = idA;
        }

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

    /** @param {number} c @param {number} ts @param {number} te @param {number} sId @param {number} eId */
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
                    if (haveGap) {
                        pushArc(c, gapTheta, evTheta[k], gapId, evId[k]);
                        haveGap = false;
                    } else {
                        seamTheta = evTheta[k];
                        seamId = evId[k];
                        haveSeam = true;
                    }
                }
                depth++;
            } else {
                depth--;
                if (depth === 0) {
                    gapTheta = evTheta[k];
                    gapId = evId[k];
                    haveGap = true;
                }
            }
        }
        if (haveGap && haveSeam) pushArc(c, gapTheta, seamTheta, gapId, seamId); // arc straddling the seam
    }

    // planar arrangement bound: the union boundary has ≤ 6·active − 12 arcs (Euler, active ≥ 3).
    // Exceeding it means the arc sweep produced spurious arcs — an internal-consistency bug.
    let active = 0;
    for (let c = 0; c < n; c++) if (!covered[c]) active++;
    if (active >= 3 && arcCount > 6 * active - 12) throw new Error('Arc count exceeds the planar 6n−12 bound.');

    return {arcCount, fullCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId, arcEndId};
}

/**
 * Stitch boundary arcs into closed rings and assemble them into the final arc topology — the `arcs()`
 * output, before any sampling. The two jobs are fused so no per-ring index bookkeeping is materialized
 * between them: each ring is built directly as its array of `[lng, lat, radius, startAngle, endAngle]` arc
 * tuples and filed into its polygon as it closes.
 *
 * **Walk.** Every CCW-oriented arc keeps the disk interior on its left, so each intersection point is the
 * *end* of exactly one arc and the *start* of exactly one other (on the other circle). We map point-ID → the
 * arc starting there, then walk end→start handoffs — an exact integer match, no geometry — until the ring
 * closes. Full-circle arcs (`startId = −1`) are standalone single-arc rings (an isolated disk's whole
 * boundary). A broken handoff (dead end or already-consumed next arc) means a degenerate input slipped past
 * dedup — an internal-consistency violation, so it throws rather than emit a corrupt ring.
 *
 * **Nest.** Each ring belongs to a polygon by *connectivity*, not geometric search: a connected component of
 * the union is one shell + zero or more holes, and a hole always belongs to its own component's shell. So
 * group rings by `scan`'s component id (shared by all arcs in a ring) — the single positive-area ring is the
 * shell, negatives are its holes — and emit one polygon per component. No point-in-ring containment test. An
 * isolated disk and an "island in a hole" are each their own component → their own polygon, exactly as RFC
 * 7946 wants (the case a naive point-in-ring nester mis-parents).
 *
 * Orientation is the sign of each ring's signed spherical area, computed analytically as we walk: geodesic
 * polygon through the arc endpoints (fan of signed spherical triangles from the first vertex) plus, per arc,
 * the segment between the small-circle arc and its chord `Δθ·(1 − cosρ) − triExcess(c, A, B)` (cap sector
 * minus geodesic triangle). A full circle's ring is the whole cap `2π·(1 − cosρ)`. Shells come out CCW
 * (positive area — already RFC 7946 winding, so sampling needs no reversal), holes CW (negative). Arc angles
 * are pre-unwrapped (`endAngle = startAngle + sweep`, `sweep ∈ (0, 2π]`); a full circle yields `[…, 0, 2π]`.
 *
 * A component missing its shell (or carrying two) is an internal-consistency violation and throws.
 *
 * @param {ReturnType<typeof build>} state
 * @param {{points: Float64Array, pointCount: number, component: Int32Array,
 *   componentCount: number}} scanResult
 * @param {{arcCount: number, arcCircle: Int32Array, arcThetaStart: Float64Array,
 *   arcThetaEnd: Float64Array, arcStartId: Int32Array, arcEndId: Int32Array}} arcResult
 * @returns {Topology}
 */
function stitch(state, scanResult, arcResult) {
    const {lng, lat, r, cx, cy, cz, cosR} = state;
    const {points, pointCount, component, componentCount} = scanResult;
    const {arcCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId, arcEndId} = arcResult;

    // Resolve the end→start handoff by point ID: the arc starting at each intersection point.
    const arcByStartId = new Int32Array(pointCount).fill(-1);
    for (let k = 0; k < arcCount; k++) if (arcStartId[k] !== -1) arcByStartId[arcStartId[k]] = k;

    // Per component, filled as rings close: the shell (positive area) and a lazily-created hole list.
    // `null` holes for the common hole-free component avoids an empty-array alloc.
    /** @type {(Ring | null)[]} */ const shellOf = new Array(componentCount).fill(null);
    /** @type {(Ring[] | null)[]} */ const holesOf = new Array(componentCount).fill(null);

    // file a just-closed ring into its component as the shell (area ≥ 0) or a hole.
    /** @param {Ring} ring @param {number} comp @param {number} area */
    const fileRing = (ring, comp, area) => {
        if (area >= 0) {
            if (shellOf[comp]) throw new Error('A connected component has more than one shell ring.');
            shellOf[comp] = ring;
        } else (holesOf[comp] ??= []).push(ring);
    };

    const visited = new Uint8Array(arcCount);
    let consumed = 0;

    for (let k0 = 0; k0 < arcCount; k0++) {
        if (visited[k0]) continue;
        const comp = component[arcCircle[k0]];

        if (arcStartId[k0] === -1) { // full circle → standalone shell ring (whole cap)
            visited[k0] = 1; consumed++;
            const c = arcCircle[k0];
            fileRing([[lng[c], lat[c], r[c], 0, TWO_PI]], comp, TWO_PI * (1 - cosR[c]));
            continue;
        }

        /** @type {Ring} */ const ring = [];
        let area = 0;
        let p0x = 0, p0y = 0, p0z = 0, havePoint0 = false;
        let k = k0;
        for (;;) {
            visited[k] = 1; consumed++;

            const sp = arcStartId[k] * 3, ep = arcEndId[k] * 3;
            const ax = points[sp], ay = points[sp + 1], az = points[sp + 2];
            const bx = points[ep], by = points[ep + 1], bz = points[ep + 2];
            const c = arcCircle[k];

            let dth = arcThetaEnd[k] - arcThetaStart[k];
            if (dth < 0) dth += TWO_PI;
            ring.push([lng[c], lat[c], r[c], arcThetaStart[k], arcThetaStart[k] + dth]);

            // segment between arc and its chord: cap sector Δθ(1−cosρ) minus geodesic triangle
            area += dth * (1 - cosR[c]);
            area -= triExcess(cx[c], cy[c], cz[c], ax, ay, az, bx, by, bz);

            // geodesic polygon through endpoints, as a fan from the ring's first vertex
            if (!havePoint0) {
                p0x = ax; p0y = ay; p0z = az; havePoint0 = true;
            } else {
                area += triExcess(p0x, p0y, p0z, ax, ay, az, bx, by, bz);
            }

            const next = arcByStartId[arcEndId[k]];
            if (next === k0) break; // ring closed
            if (next === -1 || visited[next]) {
                throw new Error('Ring failed to close — arc handoff broke (likely an undeduplicated coincident circle).');
            }
            k = next;
        }
        fileRing(ring, comp, area);
    }

    if (consumed !== arcCount) throw new Error('Not every arc was consumed exactly once while stitching rings.');

    /** @type {Topology} */ const topology = [];
    for (let comp = 0; comp < componentCount; comp++) {
        const shell = shellOf[comp];
        if (!shell) throw new Error('A connected component has no shell ring.');
        /** @type {Polygon} */ const poly = [shell];
        const holes = holesOf[comp];
        if (holes) for (const h of holes) poly.push(h);
        topology.push(poly);
    }
    return topology;
}

/**
 * Sample the arc topology into a GeoJSON `MultiPolygon`. Self-contained: it consumes only the public arc
 * shape, recomputing each arc's frame from its `[lng, lat, radius]`, so `geojson()` and an external
 * resampler walk the exact same path.
 *
 * Each arc is sampled along the exact geodesic circle `p(θ) = cosρ·c + sinρ·(cosθ·u + sinθ·v)`, then
 * projected to `[lng, lat]`. The step adapts to circle size: the largest `Δθ` whose chord stays within
 * `tolerance` km of the true arc (sagitta `r·(1−cos(Δθ/2)) ≤ tol`, so `Δθ_max = 2·acos(1 − tol/r)`), floored
 * at `minPoints` per full turn so even a tiny circle stays round. Each arc emits its start + interior
 * samples and *omits* its end (the next arc's shared start); the ring is closed by repeating its first vertex.
 *
 * Antimeridian/pole straddles are not yet split; regional data (e.g. the Ukraine workload) never triggers them.
 *
 * @param {Topology} topology
 * @param {{tolerance?: number, minPoints?: number}} [options] `tolerance`: max chord sagitta in km (default
 *   0.005 = 5 m); `minPoints`: floor on samples per full circle (default 24, keeps small circles round)
 * @returns {{type: 'MultiPolygon', coordinates: number[][][][]}}
 */
function sample(topology, options = {}) {
    const {tolerance = 0.005, minPoints = 24} = options;
    const tol = tolerance > 0 ? tolerance : 0.005;      // km, ~5 m sagitta
    const minPts = minPoints > 0 ? minPoints : 24;      // floor on points per full circle

    const coordinates = topology.map(poly => poly.map(ring => sampleRing(ring, tol, minPts)));
    return {type: 'MultiPolygon', coordinates};
}

/**
 * Sample one arc ring into a closed `[lng, lat]` vertex ring.
 * @param {Ring} ring @param {number} tol @param {number} minPts @returns {number[][]}
 */
function sampleRing(ring, tol, minPts) {
    /** @type {number[][]} */ const out = [];
    for (const [lngDeg, latDeg, rad, t0, t1] of ring) {
        // recompute the local frame at this arc's circle (mirrors `build`)
        const latR = latDeg * RAD, lngR = lngDeg * RAD;
        const cosLat = Math.cos(latR), sinLat = Math.sin(latR);
        const cosLng = Math.cos(lngR), sinLng = Math.sin(lngR);
        const cx = cosLat * cosLng, cy = cosLat * sinLng, cz = sinLat;
        const ux = -sinLng, uy = cosLng;                              // east (uz = 0)
        const vx = -sinLat * cosLng, vy = -sinLat * sinLng, vz = cosLat; // north
        const rho = rad / R, cr = Math.cos(rho), sr = Math.sin(rho);

        const dth = t1 - t0; // pre-unwrapped sweep ∈ (0, 2π]
        const maxStep = 2 * Math.acos(Math.max(-1, 1 - tol / rad));
        const segs = Math.max(1, Math.ceil(dth / maxStep), Math.ceil(minPts * dth / TWO_PI));
        const dt = dth / segs;
        for (let s = 0; s < segs; s++) {
            const th = t0 + s * dt, cth = Math.cos(th), sth = Math.sin(th);
            const px = cr * cx + sr * (cth * ux + sth * vx);
            const py = cr * cy + sr * (cth * uy + sth * vy);
            let pz = cr * cz + sr * (sth * vz); // uz = 0, so the cth·uz term drops
            if (pz > 1) pz = 1; else if (pz < -1) pz = -1; // asin domain guard
            out.push([Math.atan2(py, px) / RAD, Math.asin(pz) / RAD]);
        }
    }
    if (out.length) out.push(out[0].slice()); // close the ring (RFC 7946)
    return out;
}

/**
 * Signed area of the spherical triangle (a, b, c) of unit vectors — positive when
 * a→b→c winds CCW seen from outside. Van Oosterom–Strackee: numerically stable, no
 * `acos`, sign carried by the triple product.
 * @param {number} ax @param {number} ay @param {number} az
 * @param {number} bx @param {number} by @param {number} bz
 * @param {number} cx @param {number} cy @param {number} cz
 */
function triExcess(ax, ay, az, bx, by, bz, cx, cy, cz) {
    const crx = by * cz - bz * cy;
    const cry = bz * cx - bx * cz;
    const crz = bx * cy - by * cx;
    const triple = ax * crx + ay * cry + az * crz;
    const ab = ax * bx + ay * by + az * bz;
    const bc = bx * cx + by * cy + bz * cz;
    const ca = cx * ax + cy * ay + cz * az;
    return 2 * Math.atan2(triple, 1 + ab + bc + ca);
}

/**
 * Disjoint-set find with path halving.
 * @param {Int32Array} parent @param {number} x @returns {number} the set root of x
 */
function dsuFind(parent, x) {
    while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}
/**
 * Disjoint-set union by size.
 * @param {Int32Array} parent @param {Int32Array} setSize @param {number} a @param {number} b
 */
function dsuUnion(parent, setSize, a, b) {
    let ra = dsuFind(parent, a), rb = dsuFind(parent, b);
    if (ra === rb) return;
    if (setSize[ra] < setSize[rb]) { const t = ra; ra = rb; rb = t; }
    parent[rb] = ra; setSize[ra] += setSize[rb];
}

/** @param {Int32Array} a @param {number} cap */
function grow32(a, cap) { const g = new Int32Array(cap); g.set(a); return g; }
/** @param {Float64Array} a @param {number} cap */
function grow64(a, cap) { const g = new Float64Array(cap); g.set(a); return g; }
