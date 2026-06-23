// Fast, deterministic test suite — `npm test` (node --test).
//
// Two layers:
//   • the real OpenCelliD fixture (~23k disks) — structural invariants, a golden area
//     snapshot, and the independent point/arc oracle checks (the classify check is
//     O(covered·rimN) and too slow here, so it runs on the synthetic cases instead);
//   • hand-built synthetic cases with known topology — every check including classify,
//     run in microseconds.
//
// No Monte-Carlo, no randomness: the analytic union area is pinned to a golden constant and the
// oracle checks sweep rims at fixed angles, so a failure is a real regression, never flake.

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {build, scan, arcs, stitch, polygons} from '../index.js';
import {loadCells} from './fixtures.js';
import {checkClassify, checkPoints, checkArcs, checkArcsComplete} from './oracle.js';

const R = 6371; // mean Earth radius, km

/** Run the full pipeline and fold in the handful of scalars the tests assert on. */
function union(lng, lat, r) {
    const state = build(Float64Array.from(lng), Float64Array.from(lat), Float64Array.from(r));
    const scanResult = scan(state);
    const arcResult = arcs(state, scanResult);
    const ringResult = stitch(state, scanResult, arcResult);
    const geojson = polygons(state, scanResult, arcResult, ringResult);

    let shells = 0, holes = 0, areaSr = 0, arcsConsumed = 0;
    for (let i = 0; i < ringResult.ringCount; i++) {
        if (ringResult.ringArea[i] >= 0) shells++; else holes++;
        areaSr += ringResult.ringArea[i];
        arcsConsumed += ringResult.ringStart[i + 1] - ringResult.ringStart[i];
    }
    return {
        state, scanResult, arcResult, ringResult, geojson,
        shells, holes, arcsConsumed, areaKm2: areaSr * R * R,
        covered: scanResult.coveredCount,
        components: scanResult.componentCount,
        active: state.n - scanResult.coveredCount,
        arcCount: arcResult.arcCount,
    };
}

/** Assert a check result passes; on failure surface the first few offenders. */
function assertCheck(result) {
    assert.ok(result.pass, `${result.name}: ${result.failures.length} failure(s) — ${JSON.stringify(result.failures.slice(0, 3))}`);
}

/** Spherical-cap area of a single disk of radius `km`, in km². */
const capArea = km => 2 * Math.PI * (1 - Math.cos(km / R)) * R * R;

// The real fixture is deterministic, so compute its union once and share it across the real tests.
const cells = loadCells();
let realCache;
const realUnion = () => (realCache ??= union(cells.lng, cells.lat, cells.r));

test('real fixture — structural invariants', () => {
    const u = realUnion();
    assert.equal(u.ringResult.openRings, 0, 'every ring must close (pure-ID handoff)');
    assert.equal(u.arcsConsumed, u.arcCount, 'every arc consumed exactly once');
    assert.equal(u.shells, u.components, 'one shell per connected component');
    assert.equal(u.geojson.coordinates.length, u.components, 'one Polygon per component');
    assert.ok(u.arcCount <= Math.max(0, 6 * u.active - 12), 'planar arc-count bound');
    // current real-data topology — a coarse regression tripwire, loosen if the fixture changes
    assert.equal(u.components, 36);
    assert.equal(u.holes, 10);
});

test('real fixture — golden union area', () => {
    const u = realUnion();
    const golden = 220951.293939; // km², analytic Σ signed ring area
    assert.ok(Math.abs(u.areaKm2 - golden) / golden < 1e-6,
        `area ${u.areaKm2.toFixed(3)} km² drifted from golden ${golden} km²`);
});

test('real fixture — oracle checks (points, arcs, completeness)', () => {
    const u = realUnion();
    assertCheck(checkPoints(u.state, u.scanResult));
    assertCheck(checkArcs(u.state, u.arcResult));
    assertCheck(checkArcsComplete(u.state, u.scanResult, u.arcResult, {rimN: 64}));
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

// --- synthetic cases: known topology, all checks including classify, all instant ---

test('two overlapping circles → one shell', () => {
    const u = union([0, 0.1], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.shells, u.holes, u.covered], [1, 1, 0, 0]);
    assertCheck(checkClassify(u.state, u.scanResult));
    assertCheck(checkPoints(u.state, u.scanResult));
    assertCheck(checkArcs(u.state, u.arcResult));
    assertCheck(checkArcsComplete(u.state, u.scanResult, u.arcResult));
});

test('two disjoint circles → two components', () => {
    const u = union([0, 1], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.shells, u.holes, u.covered], [2, 2, 0, 0]);
    assert.ok(Math.abs(u.areaKm2 - 2 * capArea(10)) / u.areaKm2 < 1e-6, 'two full disks');
    assertCheck(checkArcsComplete(u.state, u.scanResult, u.arcResult));
});

test('small circle engulfed by a large one → covered', () => {
    const u = union([0, 0], [0, 0.05], [50, 1]);
    assert.deepEqual([u.components, u.shells, u.holes, u.covered], [1, 1, 0, 1]);
    assert.ok(Math.abs(u.areaKm2 - capArea(50)) / u.areaKm2 < 1e-6, 'area = the big disk');
    assertCheck(checkClassify(u.state, u.scanResult));
});

test('exact-duplicate circles are deduplicated', () => {
    const u = union([0, 0], [0, 0], [10, 10]);
    assert.deepEqual([u.components, u.shells, u.covered], [1, 1, 1]);
    assert.ok(Math.abs(u.areaKm2 - capArea(10)) / u.areaKm2 < 1e-6, 'one disk remains');
});

test('ring of circles encloses a hole', () => {
    const D = (20 / R) * 180 / Math.PI; // 20 km offset, in degrees
    const lng = [], lat = [], r = [];
    for (let k = 0; k < 8; k++) {
        const a = k * Math.PI / 4;
        lng.push(D * Math.cos(a)); lat.push(D * Math.sin(a)); r.push(9);
    }
    const u = union(lng, lat, r);
    assert.deepEqual([u.components, u.shells, u.holes], [1, 1, 1]);
    assert.equal(u.geojson.coordinates[0].length, 2, 'one Polygon with an outer ring + one hole');
    assertCheck(checkArcsComplete(u.state, u.scanResult, u.arcResult));
});
