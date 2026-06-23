import Flatbush from 'flatbush';
import {within} from 'geoflatbush';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km
const TWO_PI = 2 * Math.PI;
const COINCIDENT = 1 - 1e-13; // cos(angular dist) above which two centers are treated as co-located (§2)

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

    // Drop exact-duplicate circles up front (co-located tower rows — bit-identical lng/lat
    // are common in real data). This MUST precede the main sweep: otherwise an earlier,
    // larger circle pairs with several duplicates at once, minting bit-identical §3 points
    // under distinct IDs that desync the §5 handoff. A plain hash set keyed on the exact
    // centre does it in O(n): the first circle seen at a location wins (lowest index ⇒
    // largest radius, since circles are radius-descending) and any later one with the same
    // centre is engulfed. Only *bit-identical* centres need this — near-but-not-exact
    // coincidences mint distinct points and are caught by the §2 engulf test below, so the
    // match is exact float equality (the round() only derives hash bits, not a tolerance).
    let cap = 1; while (cap < n * 2) cap <<= 1;
    const mask = cap - 1;
    const slot = new Int32Array(cap).fill(-1);
    for (let i = 0; i < n; i++) {
        let h = (Math.imul(Math.round(lat[i] * 1e7), 0x9e3779b1) ^ Math.round(lng[i] * 1e7)) & mask;
        while (slot[h] !== -1 && (lat[slot[h]] !== lat[i] || lng[slot[h]] !== lng[i])) h = (h + 1) & mask;
        if (slot[h] === -1) slot[h] = i;                    // new location → representative
        else { covered[i] = 1; coveredCount++; }            // same centre as an earlier, larger circle
    }

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

            // proper intersection — solve p± = α·cᵢ + β·cⱼ ± γ·n for the two points.
            // n = cᵢ × cⱼ; |n|² = sin²d is computed from the cross product directly, NOT
            // as 1−cij²: for close or small circles cij≈1 and 1−cij² cancels catastrophically
            // (loses ~8 digits, throwing the points ~hundreds of m off and desyncing §4's
            // bearings between the two circles). The α/β numerators are likewise rewritten
            // as (cosRᵢ−cosRⱼ)+cosRⱼ·(1−cosd) with 1−cosd = 2sin²(d/2) to avoid the same
            // cancellation in cosRᵢ−cij·cosRⱼ. This keeps points on-circle to ~1e-15.
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

/**
 * §5. Stitch boundary arcs into closed rings, and label each ring with the data
 * §7 needs *before* sampling inflates the vertex count:
 *   - **connected-component id** (from the §3 union-find) — the component of any
 *     of the ring's arcs' circles (all arcs in a ring share one component);
 *   - **signed spherical area** — analytic closed form summed over the arcs, used
 *     for orientation (shell vs hole) under RFC 7946 winding.
 *
 * Topology: every CCW-oriented arc keeps the disk interior on its left, so each
 * §3 intersection point is the *end* of exactly one arc and the *start* of exactly
 * one other (on the other circle). We map point-ID → the arc starting there, then
 * walk end→start handoffs until the ring closes. Full-circle arcs (`startId = −1`)
 * are standalone single-arc rings (an isolated disk's whole boundary).
 *
 * Area decomposition (exact on the sphere): ring area = geodesic polygon through
 * the arc endpoints (fan of signed spherical triangles from the first vertex) plus,
 * per arc, the segment between the small-circle arc and its chord —
 * `Δθ·(1 − cosρ) − triExcess(c, A, B)` (cap sector minus geodesic triangle). A full
 * circle's ring is the whole cap, `2π·(1 − cosρ)`. Shells come out CCW (positive),
 * holes CW (negative).
 *
 * @param {State} state
 * @param {{points: Float64Array, pointCount: number, component: Int32Array}} scanResult
 * @param {{arcCount: number, arcCircle: Int32Array, arcThetaStart: Float64Array,
 *   arcThetaEnd: Float64Array, arcStartId: Int32Array, arcEndId: Int32Array}} arcResult
 * `openRings` counts rings whose end→start handoff broke; expected 0 — a nonzero count
 * signals a §2/§3/§4 inconsistency (e.g. an un-collapsed coincident circle pair).
 *
 * @returns {{ringCount: number, ringArcs: Int32Array, ringStart: Int32Array,
 *   ringComponent: Int32Array, ringArea: Float64Array, openRings: number}}
 */
export function stitch(state, scanResult, arcResult) {
    const {cx, cy, cz, cosR} = state;
    const {points, pointCount, component} = scanResult;
    const {arcCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId, arcEndId} = arcResult;

    // Resolve the end→start handoff by point ID. Each §3 intersection point is the single
    // vertex of one circle *pair*: both circles reference the *same* point ID there (§4), so
    // the arc ending at a vertex and the arc continuing from it carry identical IDs — an
    // exact integer match, no geometry. (This relies on §2 having collapsed coincident
    // duplicate circles, which would otherwise mint two IDs at one location; see scan.)
    const arcByStartId = new Int32Array(pointCount).fill(-1);
    for (let k = 0; k < arcCount; k++) if (arcStartId[k] !== -1) arcByStartId[arcStartId[k]] = k;

    const visited = new Uint8Array(arcCount);
    const ringArcs = new Int32Array(arcCount); // arc indices, grouped by ring (each arc written once)
    let w = 0;                                 // write cursor into ringArcs
    const ringStart = [0];                     // ring r spans ringArcs[ringStart[r]..ringStart[r+1])
    /** @type {number[]} */ const ringComponent = [];
    /** @type {number[]} */ const ringArea = [];
    let openRings = 0;                         // rings whose end→start handoff broke (degenerate input)

    for (let k0 = 0; k0 < arcCount; k0++) {
        if (visited[k0]) continue;
        const comp = component[arcCircle[k0]];

        if (arcStartId[k0] === -1) { // full circle → standalone ring (whole cap)
            visited[k0] = 1;
            ringArcs[w++] = k0;
            ringComponent.push(comp);
            ringArea.push(TWO_PI * (1 - cosR[arcCircle[k0]]));
            ringStart.push(w);
            continue;
        }

        let area = 0;
        let p0x = 0, p0y = 0, p0z = 0, havePoint0 = false;
        let k = k0, closed = false;
        for (;;) {
            visited[k] = 1;
            ringArcs[w++] = k;

            const sp = arcStartId[k] * 3, ep = arcEndId[k] * 3;
            const ax = points[sp], ay = points[sp + 1], az = points[sp + 2];
            const bx = points[ep], by = points[ep + 1], bz = points[ep + 2];
            const c = arcCircle[k];

            // segment between arc and its chord: cap sector Δθ(1−cosρ) minus geodesic triangle
            let dth = arcThetaEnd[k] - arcThetaStart[k];
            if (dth < 0) dth += TWO_PI;
            area += dth * (1 - cosR[c]);
            area -= triExcess(cx[c], cy[c], cz[c], ax, ay, az, bx, by, bz);

            // geodesic polygon through endpoints, as a fan from the ring's first vertex
            if (!havePoint0) { p0x = ax; p0y = ay; p0z = az; havePoint0 = true; }
            else area += triExcess(p0x, p0y, p0z, ax, ay, az, bx, by, bz);

            const next = arcByStartId[arcEndId[k]];
            if (next === k0) { closed = true; break; }
            // a dead end (−1) or an arc already consumed means the shared-ID handoff broke —
            // only happens on near-coincident circles (duplicate point IDs); stop, don't corrupt.
            if (next === -1 || visited[next]) break;
            k = next;
        }
        if (!closed) openRings++;

        ringComponent.push(comp);
        ringArea.push(area);
        ringStart.push(w);
    }

    return {
        ringCount: ringStart.length - 1,
        ringArcs,
        ringStart: Int32Array.from(ringStart),
        ringComponent: Int32Array.from(ringComponent),
        ringArea: Float64Array.from(ringArea),
        openRings,
    };
}

/**
 * §6–§7. Sample boundary arcs into vertices and assemble the GeoJSON `MultiPolygon`.
 *
 * §6 — each surviving arc is sampled along the exact geodesic circle
 *   `p(θ) = cosρ·c + sinρ·(cosθ·u + sinθ·v)` (not a chord of the centre), then projected
 *   to `[lng, lat]`. The step adapts to circle size: it is the largest `Δθ` whose chord
 *   stays within `tolerance` km of the true arc (sagitta `r·(1−cos(Δθ/2)) ≤ tol`, so
 *   `Δθ_max = 2·acos(1 − tol/r)`), floored at `minPoints` per full turn so even a tiny
 *   circle stays round (never collapses to a coarse few-gon that reads as a shrunken
 *   radius). A 30 km circle gets many points from the sagitta bound, a 100 m one gets the
 *   floor — instead of a radius-blind fixed count that over-samples the small ones.
 *   Consecutive arcs share an endpoint by §3 ID, so each arc emits its start + interior
 *   samples and *omits* its end (the next arc's start); the ring is closed by repeating v0.
 *
 * §7 — nesting is connectivity, not geometry: every ring already carries its §3
 *   component id and a signed spherical area. One `Polygon` per component — the single
 *   positive-area ring is the shell, negative-area rings are its holes. Our CCW traversal
 *   makes shells come out CCW and holes CW, which is exactly RFC 7946 winding, so no
 *   reversal is needed. An isolated disk and an island-in-a-hole are each their own
 *   component → their own `Polygon`, which is what RFC 7946 wants.
 *
 * Antimeridian/pole straddles are not yet split (Step 7 hardening); regional data
 * (e.g. the Ukraine workload) never triggers them.
 *
 * @param {State} state
 * @param {{componentCount: number}} scanResult
 * @param {{arcCircle: Int32Array, arcThetaStart: Float64Array, arcThetaEnd: Float64Array,
 *   arcStartId: Int32Array}} arcResult
 * @param {{ringCount: number, ringArcs: Int32Array, ringStart: Int32Array,
 *   ringComponent: Int32Array, ringArea: Float64Array}} ringResult
 * @param {{tolerance?: number, minPoints?: number}} [options] `tolerance`: max chord
 *   sagitta in km (default 0.005 = 5 m); `minPoints`: floor on samples per full circle
 *   (default 24, keeps small circles round)
 * @returns {{type: 'MultiPolygon', coordinates: number[][][][]}}
 */
export function polygons(state, scanResult, arcResult, ringResult, options = {}) {
    const {r: radius, cx, cy, cz, ux, uy, uz, vx, vy, vz, cosR, sinR} = state;
    const {componentCount} = scanResult;
    const {arcCircle, arcThetaStart, arcThetaEnd, arcStartId} = arcResult;
    const {ringCount, ringArcs, ringStart, ringComponent, ringArea} = ringResult;
    const tol = options.tolerance > 0 ? options.tolerance : 0.005;      // km, ~5 m sagitta
    const minPts = options.minPoints > 0 ? options.minPoints : 24;      // floor on points per full circle

    // §6 — sample each stitched ring into a closed [lng, lat] vertex ring
    const rings = new Array(ringCount);
    for (let r = 0; r < ringCount; r++) {
        const ring = [];
        for (let a = ringStart[r]; a < ringStart[r + 1]; a++) {
            const k = ringArcs[a], c = arcCircle[k];
            const cr = cosR[c], sr = sinR[c];

            // sweep angle of this arc; a full circle (startId −1) spans the whole 2π
            const full = arcStartId[k] === -1;
            const t0 = full ? 0 : arcThetaStart[k];
            let dth = full ? TWO_PI : arcThetaEnd[k] - arcThetaStart[k];
            if (dth < 0) dth += TWO_PI;

            // largest step whose chord stays within `tol` km of the arc, sized to this
            // circle's radius; emit start + interiors, skip the end (next arc's start, shared
            // by §3 ID — closure repeats v0). A full circle still needs ≥3 vertices to be valid.
            const maxStep = 2 * Math.acos(Math.max(-1, 1 - tol / radius[c]));
            // never coarser than `minPts` per full turn, so even tiny circles stay round
            const segs = Math.max(1, Math.ceil(dth / maxStep), Math.ceil(minPts * dth / TWO_PI));
            const dt = dth / segs;
            for (let s = 0; s < segs; s++) {
                const th = t0 + s * dt, cth = Math.cos(th), sth = Math.sin(th);
                const px = cr * cx[c] + sr * (cth * ux[c] + sth * vx[c]);
                const py = cr * cy[c] + sr * (cth * uy[c] + sth * vy[c]);
                let pz = cr * cz[c] + sr * (cth * uz[c] + sth * vz[c]);
                if (pz > 1) pz = 1; else if (pz < -1) pz = -1; // asin domain guard
                ring.push([Math.atan2(py, px) / RAD, Math.asin(pz) / RAD]);
            }
        }
        if (ring.length) ring.push(ring[0].slice()); // close the ring (RFC 7946)
        rings[r] = ring;
    }

    // §7 — group rings by component: positive-area ring is the shell, negatives are its holes
    const shellOf = new Int32Array(componentCount).fill(-1);
    /** @type {number[][]} */ const holesOf = [];
    for (let comp = 0; comp < componentCount; comp++) holesOf.push([]);
    for (let r = 0; r < ringCount; r++) {
        const comp = ringComponent[r];
        if (ringArea[r] >= 0) shellOf[comp] = r; else holesOf[comp].push(r);
    }

    const coordinates = [];
    for (let comp = 0; comp < componentCount; comp++) {
        const sr = shellOf[comp];
        if (sr === -1) continue; // a component with no shell would be a §5 inconsistency
        const poly = [rings[sr]];
        for (const hr of holesOf[comp]) poly.push(rings[hr]);
        coordinates.push(poly);
    }

    return {type: 'MultiPolygon', coordinates};
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

/** @param {Int32Array} a @param {number} cap */
function grow32(a, cap) { const g = new Int32Array(cap); g.set(a); return g; }
/** @param {Float64Array} a @param {number} cap */
function grow64(a, cap) { const g = new Float64Array(cap); g.set(a); return g; }
