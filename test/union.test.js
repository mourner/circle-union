// Fast, deterministic test suite — `npm test` (node --test).
//
// End-to-end and black-box: every assertion is made against the PUBLIC CircleUnion API —
// `arcs()` (exact topology) and `finish()` (sampled GeoJSON) — never against internal
// pipeline state. Two layers:
//   • the real OpenCelliD fixture (~23k disks) — topology shape, a golden area snapshot,
//     GeoJSON well-formedness, and the independent membership-oracle check;
//   • hand-built synthetic cases with known topology, run in microseconds.
//
// Cheap structural invariants (every ring closes, every arc consumed once, one shell per
// component, arc count ≤ 6n−12) are no longer asserted here — they are runtime throws
// inside the pipeline, so any violation makes `arcs()`/`finish()` throw.
//
// No Monte-Carlo, no randomness: the area is pinned to a golden constant and the oracle
// sweeps rims at fixed angles, so a failure is a real regression, never flake.

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {CircleUnion} from '../index.js';
import {loadCells} from './fixtures.js';
import {checkTopology} from './oracle.js';

const RAD = Math.PI / 180;
const R = 6371; // mean Earth radius, km

/** Build a union from plain coordinate arrays and return {input, u, arcs, geojson, …}. */
function union(lng, lat, r, options) {
    const u = new CircleUnion(lng.length);
    for (let i = 0; i < lng.length; i++) u.add(lng[i], lat[i], r[i]);
    const arcs = u.arcs();
    const geojson = u.finish(options);
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
    const golden = 220920.667024; // km², sampled finish() output at default tolerance
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
    assert.deepEqual(u.finish(), {type: 'MultiPolygon', coordinates: []});
});

test('add past the reserved count throws', () => {
    const u = new CircleUnion(1);
    u.add(0, 0, 10);
    assert.throws(() => u.add(1, 1, 10), /reserved/);
});
