// Benchmark the full union pipeline on the real OpenCelliD fixture — `npm run bench`.
//
// Times the common public path end to end: construct a CircleUnion, add every circle, then
// `finish()` to GeoJSON (build + arc topology + sampling — sampling is a tiny tail, so this
// also stands in for `arcs()`). Warmup runs let the JIT settle and warm caches; we report the
// best (min) measured time — the most stable estimate of true cost on a noisy machine — plus
// the median. To break the total down during development, drop a temporary `console.time`
// around a stage; we deliberately keep no per-stage exports so the public surface stays just
// `CircleUnion`.
import {CircleUnion} from '../index.js';
import {loadCells} from '../test/fixtures.js';

const WARMUP = 5;
const RUNS = 30;

const {n, lng, lat, r} = loadCells();

const times = [];
for (let i = 0; i < WARMUP + RUNS; i++) {
    const t0 = performance.now();
    const u = new CircleUnion(n);
    for (let c = 0; c < n; c++) u.add(lng[c], lat[c], r[c]);
    u.finish();
    if (i >= WARMUP) times.push(performance.now() - t0);
}
times.sort((a, b) => a - b);

console.log(`circle-union benchmark — ${n.toLocaleString()} circles, ${WARMUP} warmup + ${RUNS} runs\n`);
console.log(`finish() → GeoJSON: ${times[0].toFixed(2)}ms min, ${times[times.length >> 1].toFixed(2)}ms median`);
