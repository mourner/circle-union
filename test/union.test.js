// Fast, deterministic test suite — `npm test` (node --test).
//
// End-to-end and black-box: every assertion is made against the PUBLIC CircleUnion API — `arcs()` (exact
// topology) and `geojson()` (sampled GeoJSON) — never against internal pipeline state. Two layers:
//   • the real OpenCelliD fixture (~23k disks) — topology shape, a golden area snapshot, GeoJSON
//     well-formedness, and the independent membership-oracle check;
//   • hand-built synthetic cases with known topology, run in microseconds.
//
// Cheap structural invariants (every ring closes, every arc consumed once, one shell per component, arc
// count ≤ 6n−12) are no longer asserted here — they are runtime throws inside the pipeline, so any violation
// makes `arcs()`/`geojson()` throw.
//
// No Monte-Carlo, no randomness: the area is pinned to a golden constant and the oracle sweeps rims at fixed
// angles, so a failure is a real regression, never flake.

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import Flatbush from 'flatbush';
import {within} from 'geoflatbush';
import {CircleUnion} from '../index.js';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km
const TWO_PI = 2 * Math.PI;

// The OpenCelliD Ukraine cell-tower sample (~23k disks): the workload the library was built for. Trimmed to
// the only columns we use — `lon,lat,range_m` — one tower per line.
function loadCells() {
    const lines = readFileSync(new URL('./fixtures/ukraine-cell-id.csv', import.meta.url), 'utf8').split('\n');
    const lng = new Float64Array(lines.length), lat = new Float64Array(lines.length), r = new Float64Array(lines.length);
    let n = 0;
    for (const line of lines) {
        if (!line) continue;
        const c = line.split(',');
        lng[n] = +c[0]; lat[n] = +c[1]; r[n] = +c[2] / 1000; // m → km
        n++;
    }
    return {lng: lng.subarray(0, n), lat: lat.subarray(0, n), r: r.subarray(0, n)};
}

/** Build a union from plain coordinate arrays and return {input, u, arcs, geojson, …}. */
function union(lng, lat, r, options) {
    const u = new CircleUnion(lng.length);
    for (let i = 0; i < lng.length; i++) u.add(lng[i], lat[i], r[i]);
    const arcs = u.arcs();
    const geojson = u.geojson(options);
    let holes = 0;
    for (const poly of arcs) holes += poly.length - 1;
    return {
        input: {lng, lat, r}, u, arcs, geojson,
        components: arcs.length, // one polygon (shell + holes) per connected component
        holes,
        areaKm2: geojsonAreaKm2(geojson),
    };
}

/** Signed area (steradians) of a closed [lng,lat] ring — standard spherical formula. */
function ringSteradians(ring) {
    let total = 0;
    for (let i = 0, len = ring.length - 1; i < len; i++) {
        const [lng1, lat1] = ring[i], [lng2, lat2] = ring[i + 1];
        total += (lng2 - lng1) * RAD * (2 + Math.sin(lat1 * RAD) + Math.sin(lat2 * RAD));
    }
    return total / 2;
}

/** Union area (km²) of a sampled GeoJSON MultiPolygon: per polygon, shell minus its holes. */
function geojsonAreaKm2(geojson) {
    let sr = 0;
    for (const poly of geojson.coordinates) {
        sr += Math.abs(ringSteradians(poly[0]));
        for (let h = 1; h < poly.length; h++) sr -= Math.abs(ringSteradians(poly[h]));
    }
    return sr * R * R;
}

/** Assert an oracle check result passes; on failure surface the first few offenders. */
function assertCheck(result) {
    assert.ok(result.pass, `${result.name}: ${result.failures.length} failure(s) — ${JSON.stringify(result.failures.slice(0, 3))}`);
}

/** Spherical-cap area of a single disk of radius `km`, in km². */
const capArea = km => 2 * Math.PI * (1 - Math.cos(km / R)) * R * R;

// The real fixture is deterministic, so compute its union once and share it across the real tests.
const cells = loadCells();
let realCache;
const realUnion = () => (realCache ??= union(cells.lng, cells.lat, cells.r));

test('real fixture — topology shape', () => {
    const u = realUnion();
    // current real-data topology — a coarse regression tripwire, loosen if the fixture changes
    assert.equal(u.components, 36);
    assert.equal(u.holes, 10);
    // arc shape: every arc is [lng, lat, r, θ0, θ1] with a positive, ≤2π pre-unwrapped sweep
    for (const poly of u.arcs) for (const ring of poly) for (const arc of ring) {
        assert.equal(arc.length, 5);
        const sweep = arc[4] - arc[3];
        assert.ok(sweep > 0 && sweep <= 2 * Math.PI + 1e-9, 'sweep in (0, 2π]');
    }
});

test('real fixture — golden union area', () => {
    const u = realUnion();
    const golden = 220920.667024; // km², sampled geojson() output at default tolerance
    assert.ok(Math.abs(u.areaKm2 - golden) / golden < 1e-6,
        `area ${u.areaKm2.toFixed(6)} km² drifted from golden ${golden} km²`);
});

test('real fixture — membership oracle', () => {
    const u = realUnion();
    assertCheck(checkTopology(u.input, u.arcs, {rimN: 64}));
});

test('real fixture — GeoJSON is well-formed', () => {
    const u = realUnion();
    assert.equal(u.geojson.type, 'MultiPolygon');
    for (const poly of u.geojson.coordinates) {
        assert.ok(poly.length >= 1, 'a Polygon has at least an outer ring');
        for (const ring of poly) {
            assert.ok(ring.length >= 4, 'a ring has ≥4 positions (closed triangle minimum)');
            assert.deepEqual(ring[0], ring[ring.length - 1], 'ring is closed');
            for (const [lng, lat] of ring) {
                assert.ok(lng >= -180 && lng <= 180 && lat >= -90 && lat <= 90, 'position in range');
            }
        }
    }
});

// --- synthetic cases: known topology, oracle on every one, all instant ---

test('two overlapping circles → one shell', () => {
    const u = union([0, 0.1], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.holes], [1, 0]);
    assertCheck(checkTopology(u.input, u.arcs));
});

test('two disjoint circles → two components', () => {
    const u = union([0, 1], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.holes], [2, 0]);
    assert.ok(Math.abs(u.areaKm2 - 2 * capArea(10)) / u.areaKm2 < 2e-3, 'two full disks');
    assertCheck(checkTopology(u.input, u.arcs));
});

test('small circle engulfed by a large one → covered', () => {
    const u = union([0, 0], [0, 0.05], [50, 1]);
    assert.deepEqual([u.components, u.holes], [1, 0]);
    assert.ok(Math.abs(u.areaKm2 - capArea(50)) / u.areaKm2 < 2e-3, 'area = the big disk');
    assertCheck(checkTopology(u.input, u.arcs));
});

test('exact-duplicate circles are deduplicated', () => {
    const u = union([0, 0], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.holes], [1, 0]);
    assert.ok(Math.abs(u.areaKm2 - capArea(10)) / u.areaKm2 < 2e-3, 'one disk remains');
    assertCheck(checkTopology(u.input, u.arcs));
});

test('ring of circles encloses a hole', () => {
    const D = (20 / R) * 180 / Math.PI; // 20 km offset, in degrees
    const lng = [], lat = [], r = [];
    for (let k = 0; k < 8; k++) {
        const a = k * Math.PI / 4;
        lng.push(D * Math.cos(a)); lat.push(D * Math.sin(a)); r.push(9);
    }
    const u = union(lng, lat, r);
    assert.deepEqual([u.components, u.holes], [1, 1]);
    assert.equal(u.arcs[0].length, 2, 'one Polygon with an outer ring + one hole');
    assertCheck(checkTopology(u.input, u.arcs));
});

test('empty union → empty output', () => {
    const u = new CircleUnion(0);
    assert.deepEqual(u.arcs(), []);
    assert.deepEqual(u.geojson(), {type: 'MultiPolygon', coordinates: []});
});

test('add past the reserved count throws', () => {
    const u = new CircleUnion(1);
    u.add(0, 0, 10);
    assert.throws(() => u.add(1, 1, 10), /reserved/);
});

// --- independent membership oracle -----------------------------------------------------------
//
// The whole algorithm is validated against ONE independent ground truth: a brute-force membership oracle. A
// point P (unit vector) is in the union iff it lies inside some disk, i.e. dot(P, centerᵢ) ≥ cosRᵢ for some
// i — an O(n) test sharing no code with the arc algorithm, so a bug in one is very unlikely to be mirrored
// in the other.
//
// BLACK-BOX BY CONSTRUCTION. `checkTopology` is handed the test's *own input circles* and the *public
// `arcs()` topology* — nothing internal. It rebuilds everything from the input (its own Flatbush, its own
// per-circle frames), then judges the public output against it: for each input circle, sweep its rim and
// classify every sample with the oracle (excluding the circle and any exact duplicates). The exposed runs —
// points on the union boundary — must coincide *exactly* with that circle's emitted arcs (grouped by its
// (lng, lat, r) key). This catches a covered circle that kept an arc, a spurious arc, and a *missing* arc
// (the precision corner the per-arc check is blind to); since arc endpoints are the intersection points, a
// tight-ε match also validates the solve. Deterministic (fixed-angle rim sweep): a failure is a real
// regression, never flake.

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
function checkTopology(input, topology, {rimN = 256, eps = 1e-9} = {}) {
    const {lng, lat, r} = input;
    const n = lng.length;

    // Independent per-circle geometry, recomputed from the raw input (mirrors `build` but shares no code
    // with it): center unit vector, east/north frame, cos/sin of angular radius.
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
        // Disks that can bury any of c's rim have centers within r_c + maxR — gather once with the index,
        // then keep only those that actually *reach* c's rim: disk i can cover a rim point of c only if the
        // two disks touch, dist(centerᵢ, centerc) ≤ ρ_c + ρ_i, i.e. dot(centers) ≥ cos(ρ_c+ρ_i) =
        // cosR_c·cosR_i − sinR_c·sinR_i. This prunes the lone large outlier and far small disks, so the
        // per-sample loop stays short.
        const cand = within(index, lng[c], lat[c], r[c] + maxR + 1e-6);
        const reach = [];
        for (let q = 0; q < cand.length; q++) {
            const i = cand[q];
            // self + exact duplicates (numeric compare — building a string key here, once per candidate
            // across every circle, dominated the whole check)
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
