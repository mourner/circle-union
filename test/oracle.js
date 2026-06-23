// Correctness net for CircleUnion (dev/test only — not shipped).
//
// The whole algorithm is validated against ONE independent ground truth: a brute-force
// membership oracle. A point P (unit vector) is in the union iff it lies inside some disk,
// i.e. dot(P, centerᵢ) ≥ cosRᵢ for some i — an O(n) loop sharing no code with the arc
// algorithm, so a bug in one is very unlikely to be mirrored in the other.
//
// BLACK-BOX BY CONSTRUCTION. This module is handed the test's *own input circles* and the
// *public `arcs()` topology* — nothing internal. It rebuilds everything it needs from the
// input (its own Flatbush, its own per-circle frames), recomputing the geometry from
// scratch, then judges the public output against it.
//
// A single per-circle check subsumes the old classify / soundness / completeness suite:
// for each input circle, sweep its rim and classify every sample with the oracle (excluding
// the circle and any exact duplicates of it). The exposed runs — points on the union
// boundary — must coincide *exactly* with that circle's emitted arcs (grouped by its
// (lng, lat, r) key). This catches a covered circle that kept an arc, a spurious arc, and a
// *missing* arc (the §3-precision corner the per-arc check is blind to). Because arc
// endpoints are the §3 intersection points, a tight-ε match also validates the solve.
//
// The check is DETERMINISTIC (fixed-angle rim sweep, no RNG): a failure is a real
// regression, never flake. It returns {name, total, failures, pass}; failures carry enough
// context to locate the offending circle. Nothing here throws — the caller decides.

import Flatbush from 'flatbush';
import {within} from 'geoflatbush';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km
const TWO_PI = 2 * Math.PI;

/** Stable key grouping a circle and its exact duplicates. */
const key = (lng, lat, r) => `${lng},${lat},${r}`;

/** Does angle θ (mod 2π) fall inside any pre-unwrapped [a0, a1] interval (sweep ≤ 2π)? */
function inAnyArc(intervals, theta) {
    for (const [a0, a1] of intervals) {
        let t = theta;
        while (t < a0) t += TWO_PI;
        while (t >= a0 + TWO_PI) t -= TWO_PI;
        if (t <= a1) return true;
    }
    return false;
}

/**
 * Validate the public arc topology against the independent membership oracle.
 *
 * @param {{lng: ArrayLike<number>, lat: ArrayLike<number>, r: ArrayLike<number>}} input
 *   the test's own circles (any order; duplicates allowed)
 * @param {Array<Array<Array<number[]>>>} topology  what `CircleUnion.arcs()` returned
 * @param {{rimN?: number, eps?: number}} [opts] `rimN` rim samples per circle (default 256);
 *   `eps` on-boundary band in dot-margin units (default 1e-9)
 * @returns {{name: string, total: number, failures: object[], pass: boolean}}
 */
export function checkTopology(input, topology, {rimN = 256, eps = 1e-9} = {}) {
    const {lng, lat, r} = input;
    const n = lng.length;

    // Independent per-circle geometry, recomputed from the raw input (mirrors §0 setup but
    // shares no code with it): center unit vector, east/north frame, cos/sin of angular radius.
    const cx = new Float64Array(n), cy = new Float64Array(n), cz = new Float64Array(n);
    const ux = new Float64Array(n), uy = new Float64Array(n);                       // east (uz = 0)
    const vx = new Float64Array(n), vy = new Float64Array(n), vz = new Float64Array(n); // north
    const cosR = new Float64Array(n), sinR = new Float64Array(n);
    const index = new Flatbush(n);
    let maxR = 0;
    for (let i = 0; i < n; i++) {
        const latR = lat[i] * RAD, lngR = lng[i] * RAD;
        const cosLat = Math.cos(latR), sinLat = Math.sin(latR);
        const cosLng = Math.cos(lngR), sinLng = Math.sin(lngR);
        cx[i] = cosLat * cosLng; cy[i] = cosLat * sinLng; cz[i] = sinLat;
        ux[i] = -sinLng; uy[i] = cosLng;
        vx[i] = -sinLat * cosLng; vy[i] = -sinLat * sinLng; vz[i] = cosLat;
        const rho = r[i] / R;
        cosR[i] = Math.cos(rho); sinR[i] = Math.sin(rho);
        if (r[i] > maxR) maxR = r[i];
        index.add(lng[i], lat[i]);
    }
    index.finish();

    // Group the public arcs by circle key → list of [startAngle, endAngle] intervals.
    /** @type {Map<string, Array<[number, number]>>} */ const arcsByKey = new Map();
    for (const poly of topology) for (const ring of poly) for (const [alng, alat, ar, a0, a1] of ring) {
        const k = key(alng, alat, ar);
        let list = arcsByKey.get(k);
        if (!list) arcsByKey.set(k, list = []);
        list.push([a0, a1]);
    }

    const p = [0, 0, 0];
    const failures = [];
    let total = 0;
    const seen = new Set();
    for (let c = 0; c < n; c++) {
        const k = key(lng[c], lat[c], r[c]);
        if (seen.has(k)) continue; // a duplicate shares its representative's rim and arcs
        seen.add(k);
        total++;

        const intervals = arcsByKey.get(k) || [];
        // Disks that can bury any of c's rim have centers within r_c + maxR — gather once with
        // the index, then keep only those that actually *reach* c's rim: disk i can cover a rim
        // point of c only if the two disks touch, dist(centerᵢ, centerc) ≤ ρ_c + ρ_i, i.e.
        // dot(centers) ≥ cos(ρ_c+ρ_i) = cosR_c·cosR_i − sinR_c·sinR_i. This prunes the lone large
        // outlier and far small disks, so the per-sample loop stays short.
        const cand = within(index, lng[c], lat[c], r[c] + maxR + 1e-6);
        const reach = [];
        for (let q = 0; q < cand.length; q++) {
            const i = cand[q];
            // self + exact duplicates (numeric compare — building a string key here, once per
            // candidate across every circle, dominated the whole check)
            if (lng[i] === lng[c] && lat[i] === lat[c] && r[i] === r[c]) continue;
            const dotc = cx[c] * cx[i] + cy[c] * cy[i] + cz[c] * cz[i];
            if (dotc >= cosR[c] * cosR[i] - sinR[c] * sinR[i] - 1e-12) reach.push(i);
        }

        let missing = 0, spurious = 0, firstTheta = 0;
        for (let s = 0; s < rimN; s++) {
            const theta = (s / rimN) * TWO_PI;
            const ct = Math.cos(theta), st = Math.sin(theta);
            p[0] = cosR[c] * cx[c] + sinR[c] * (ct * ux[c] + st * vx[c]);
            p[1] = cosR[c] * cy[c] + sinR[c] * (ct * uy[c] + st * vy[c]);
            p[2] = cosR[c] * cz[c] + sinR[c] * (st * vz[c]); // uz = 0

            // signed interior margin against every reaching disk
            let m = -Infinity;
            for (let q = 0; q < reach.length; q++) {
                const i = reach[q];
                const d = p[0] * cx[i] + p[1] * cy[i] + p[2] * cz[i] - cosR[i];
                if (d > m) m = d;
            }

            if (m > -eps && m < eps) continue;       // on a vertex — ambiguous, skip
            const exposed = m < 0;                    // outside all others → on the union boundary
            const inArc = inAnyArc(intervals, theta);
            if (exposed && !inArc) {
                if (!missing && !spurious) firstTheta = theta;
                missing++;
            } else if (!exposed && inArc) {
                if (!missing && !spurious) firstTheta = theta;
                spurious++;
            }
        }
        if (missing || spurious) failures.push({circle: c, key: k, missing, spurious, firstTheta});
    }

    return {name: 'topology', total, failures, pass: failures.length === 0};
}
