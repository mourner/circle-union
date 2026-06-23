// Correctness net for the union-of-disks pipeline (dev/test only — not shipped).
//
// The whole algorithm is validated against ONE independent ground truth: a brute-force
// membership oracle. A point P (unit vector) is in the union iff it lies inside some disk,
// i.e. dot(P, centerᵢ) ≥ cosRᵢ for some i — an O(n) loop sharing no code with the arc
// algorithm, so a bug in one is very unlikely to be mirrored in the other. Every pipeline
// stage emits an invariant checkable against that oracle, and all checks here are
// DETERMINISTIC (fixed-angle rim sweeps, no RNG):
//
//   classify  — every covered circle's rim lies inside the union of the others.
//   points    — each intersection point lies on both parent circles' rims (|P| = 1).
//   arcs      — (soundness) each arc's midpoint is on the boundary, buried by no other disk.
//   complete  — (completeness) sweeping each active circle's rim, the exposed runs coincide
//               exactly with its emitted arcs — catches a *missing* arc the per-arc check can't.
//
// (Area is validated separately, in the test, as a golden snapshot of the analytic Σ ringArea —
// independent random Monte-Carlo estimation lived here once but was dropped for determinism.)
//
// Each check returns {name, total, failures, pass}; failures carry enough context to locate
// the offending circle/point/arc. Nothing here throws — callers decide how to report.

import {within} from 'geoflatbush';

const RAD = Math.PI / 180;

/**
 * Build the membership oracle: a closure `margin(x,y,z, exclude?)` giving the signed interior
 * margin of a point against the union — max over disks of (dot − cosR). > 0 strictly inside some
 * disk · ≈ 0 on the boundary · < 0 outside. Candidate disks are narrowed with the flatbush index
 * (centers within maxRadius of the query); the exact dot/cosR test still decides membership, so
 * the answer is identical to the O(n) brute force, just without the n² blowup. The index is the
 * one in `state` — independent of the arc-stitching code under test.
 * @param {any} state
 * @returns {(x:number, y:number, z:number, exclude?:Set<number>) => number}
 */
export function makeOracle(state) {
    const {n, r, cx, cy, cz, cosR, index} = state;
    let maxR = 0;
    for (let i = 0; i < n; i++) if (r[i] > maxR) maxR = r[i];
    // tiny pad so a point exactly on a rim isn't missed by floating-point in the index query
    const queryR = maxR * (1 + 1e-9) + 1e-6;
    return function margin(x, y, z, exclude) {
        const plng = Math.atan2(y, x) / RAD;
        const plat = Math.asin(z < -1 ? -1 : z > 1 ? 1 : z) / RAD;
        const cand = within(index, plng, plat, queryR);
        let best = -Infinity;
        for (let c = 0; c < cand.length; c++) {
            const i = cand[c];
            if (exclude && exclude.has(i)) continue;
            const m = x * cx[i] + y * cy[i] + z * cz[i] - cosR[i];
            if (m > best) best = m;
        }
        return best;
    };
}

/** Point on circle `c` at angle θ, offset radially by `drho` radians (drho=0 → on the rim). */
function rimPoint(state, c, theta, drho, out) {
    const {cx, cy, cz, ux, uy, uz, vx, vy, vz, cosR, sinR} = state;
    const rho = Math.atan2(sinR[c], cosR[c]) + drho;
    const k = Math.cos(rho), s = Math.sin(rho);
    const ct = Math.cos(theta), st = Math.sin(theta);
    const tx = ux[c] * ct + vx[c] * st;
    const ty = uy[c] * ct + vy[c] * st;
    const tz = uz[c] * ct + vz[c] * st;
    out[0] = cx[c] * k + tx * s;
    out[1] = cy[c] * k + ty * s;
    out[2] = cz[c] * k + tz * s;
}

/**
 * Classification soundness. Every circle flagged `covered` must be globally redundant: its entire rim lies
 * inside the union of the OTHER disks, so dropping it cannot change the union. `rimN` sample points
 * per circle; `eps` is the on-boundary band. (The converse — that an *active* circle is exposed —
 * is NOT asserted: a circle can be jointly buried by several neighbours yet engulfed by none, so it
 * stays active but contributes no arc. That case is legitimate, and the area check would catch any
 * real coverage gap globally.)
 */
export function checkClassify(state, scanResult, {rimN = 256, eps = 1e-9} = {}) {
    const {n, lng, lat, r, cx, cy, cz, cosR, index} = state;
    const {covered} = scanResult;
    let maxR = 0;
    for (let i = 0; i < n; i++) if (r[i] > maxR) maxR = r[i];
    const p = [0, 0, 0];
    const failures = [];
    let total = 0;
    for (let c = 0; c < n; c++) {
        if (!covered[c]) continue;
        total++;
        // disks that can cover any of c's rim have centers within r_c + maxR — gather them once,
        // then test all rim points against this small local list (no per-point index query).
        const cand = within(index, lng[c], lat[c], r[c] + maxR + 1e-6);
        let worst = Infinity, worstTheta = 0;
        for (let s = 0; s < rimN; s++) {
            const th = (s / rimN) * 2 * Math.PI;
            rimPoint(state, c, th, 0, p);
            let m = -Infinity;
            for (let q = 0; q < cand.length; q++) {
                const i = cand[q];
                if (i === c) continue;
                const d = p[0] * cx[i] + p[1] * cy[i] + p[2] * cz[i] - cosR[i];
                if (d > m) m = d;
            }
            if (m < worst) { worst = m; worstTheta = th; }
        }
        if (worst < -eps) failures.push({circle: c, kind: 'covered-but-exposed', margin: worst, theta: worstTheta});
    }
    return {name: 'classify', total, failures, pass: failures.length === 0};
}

/**
 * Geometry of the intersection solve. Every intersection point must be a true unit-sphere point
 * lying on BOTH parent circles' rims: |P| = 1, dot(P, centerᵢ) = cosRᵢ and dot(P, centerⱼ) = cosRⱼ.
 * This checks the plane∩plane∩sphere solve in isolation. The solve deliberately keeps both roots of
 * each pair, including the one buried inside a third disk — selecting which roots bound surviving
 * arcs is `arcs`'s job (see checkArcs), not a property of the points themselves.
 */
export function checkPoints(state, scanResult, {eps = 1e-9} = {}) {
    const {cx, cy, cz, cosR} = state;
    const {points, pointCount, pairs, pairCount} = scanResult;
    // map each point ID → its two parent circles
    const parentA = new Int32Array(pointCount).fill(-1);
    const parentB = new Int32Array(pointCount).fill(-1);
    for (let pi = 0; pi < pairCount; pi++) {
        const i = pairs[pi * 3], j = pairs[pi * 3 + 1], base = pairs[pi * 3 + 2];
        parentA[base] = i; parentB[base] = j;
        parentA[base + 1] = i; parentB[base + 1] = j;
    }
    const failures = [];
    for (let id = 0; id < pointCount; id++) {
        const x = points[id * 3], y = points[id * 3 + 1], z = points[id * 3 + 2];
        const a = parentA[id], b = parentB[id];
        const unit = Math.abs(Math.hypot(x, y, z) - 1);
        const dA = Math.abs(x * cx[a] + y * cy[a] + z * cz[a] - cosR[a]);
        const dB = Math.abs(x * cx[b] + y * cy[b] + z * cz[b] - cosR[b]);
        const err = Math.max(unit, dA, dB);
        if (err > eps) failures.push({point: id, parents: [a, b], unit, onA: dA, onB: dB});
    }
    return {name: 'points', total: pointCount, failures, pass: failures.length === 0};
}

/**
 * Arc soundness. Every emitted boundary arc is genuinely on the union boundary, not buried:
 * its midpoint is interior to no disk other than its own circle (margin excl self ≤ eps).
 * Its complement — no exposed rim left without an arc (completeness) — is checked authoritatively
 * by `checkArcsComplete`'s rim sweep.
 */
export function checkArcs(state, arcResult, {eps = 1e-9} = {}) {
    const {arcCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId} = arcResult;
    const oracle = makeOracle(state);
    const TWO_PI = 2 * Math.PI;
    const p = [0, 0, 0];
    const failures = [];
    for (let k = 0; k < arcCount; k++) {
        const c = arcCircle[k];
        const self = new Set([c]);
        let dth = arcThetaEnd[k] - arcThetaStart[k];
        if (arcStartId[k] === -1) dth = TWO_PI;       // full circle
        else if (dth < 0) dth += TWO_PI;
        const mid = arcStartId[k] === -1 ? 0 : arcThetaStart[k] + dth / 2;

        rimPoint(state, c, mid, 0, p);
        const mIn = oracle(p[0], p[1], p[2], self);
        if (mIn > eps) failures.push({arc: k, circle: c, kind: 'midpoint-buried', margin: mIn});
    }
    return {name: 'arcs', total: arcCount, failures, pass: failures.length === 0};
}

/**
 * Arc completeness. The per-arc soundness check above can only judge arcs that *exist* — it is
 * structurally blind to a boundary segment `arcs` failed to emit. This check catches the missing arc:
 * for every active (non-covered) circle, sweep its rim and classify each sample with the oracle
 * (excluding self) — `margin < −eps` ⇒ exposed (on the union boundary), `> +eps` ⇒ buried, else
 * on a vertex (skipped). The exposed samples must coincide *exactly* with the circle's emitted
 * arcs: an exposed sample inside no arc is a `missing-arc` (the 3-way near-coincidence corner), a
 * buried sample inside an arc is a `spurious-arc`. Reports one failure per inconsistent circle.
 */
export function checkArcsComplete(state, scanResult, arcResult, {rimN = 256, eps = 1e-9} = {}) {
    const {n} = state;
    const {covered} = scanResult;
    const {arcCount, arcCircle, arcThetaStart, arcThetaEnd, arcStartId} = arcResult;
    const oracle = makeOracle(state);
    const TWO_PI = 2 * Math.PI;

    // gather emitted arcs per circle; `full` marks a whole-circle arc (always inside)
    /** @type {Array<Array<[number, number]>>} */ const ivals = Array.from({length: n}, () => []);
    const full = new Uint8Array(n);
    for (let k = 0; k < arcCount; k++) {
        const c = arcCircle[k];
        if (arcStartId[k] === -1) { full[c] = 1; continue; }
        ivals[c].push([arcThetaStart[k], arcThetaEnd[k]]);
    }
    const inAnyArc = (c, theta) => {
        if (full[c]) return true;
        for (const [s, e] of ivals[c]) {
            let d = (theta - s) % TWO_PI; if (d < 0) d += TWO_PI;
            let len = (e - s) % TWO_PI; if (len < 0) len += TWO_PI;
            if (d <= len) return true;
        }
        return false;
    };

    const p = [0, 0, 0];
    const failures = [];
    let total = 0;
    for (let c = 0; c < n; c++) {
        if (covered[c]) continue;
        total++;
        const self = new Set([c]);
        let missing = 0, spurious = 0, firstTheta = 0;
        for (let s = 0; s < rimN; s++) {
            const theta = (s / rimN) * TWO_PI;
            rimPoint(state, c, theta, 0, p);
            const m = oracle(p[0], p[1], p[2], self);
            if (m > -eps && m < eps) continue;          // on a vertex — ambiguous, skip
            const exposed = m < 0;
            const inArc = inAnyArc(c, theta);
            if (exposed && !inArc) {
                if (!missing && !spurious) firstTheta = theta; missing++;
            } else if (!exposed && inArc) {
                if (!missing && !spurious) firstTheta = theta; spurious++;
            }
        }
        if (missing || spurious) failures.push({circle: c, missing, spurious, firstTheta});
    }
    return {name: 'complete', total, failures, pass: failures.length === 0};
}

/** Run all deterministic stage checks and return them as an array; `format` renders one line each. */
export function runAll(state, scanResult, arcResult, opts = {}) {
    return [
        checkClassify(state, scanResult, opts.classify),
        checkPoints(state, scanResult, opts.points),
        checkArcs(state, arcResult, opts.arcs),
        checkArcsComplete(state, scanResult, arcResult, opts.arcsComplete),
    ];
}

/** Format one check result as a single status line. */
export function format(result) {
    const tag = result.pass ? 'OK  ' : 'FAIL';
    let line = `  [${tag}] ${result.name.padEnd(12)} ${result.total - result.failures.length}/${result.total} ok`;
    for (const f of result.failures.slice(0, 5)) line += `\n           ↳ ${JSON.stringify(f)}`;
    if (result.failures.length > 5) line += `\n           ↳ … ${result.failures.length - 5} more`;
    return line;
}
