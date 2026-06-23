// Benchmark the pipeline on the real OpenCelliD fixture — `npm run bench`.
//
// Each stage is timed independently with warmup runs (to let the JIT settle and warm caches)
// followed by measured runs; we report the best (min) time per stage, which is the most stable
// estimate of true cost on a noisy machine, plus the median. `build` is re-run every iteration
// because the later stages consume its allocations.
import {_stages} from '../index.js';
import {loadCells} from '../test/fixtures.js';

const {build, scan, arcs, stitch, assemble, sample} = _stages;

const WARMUP = 5;
const RUNS = 30;

const {n, lng, lat, r} = loadCells();

// Each stage takes the prior stage's output; the runner rebuilds the inputs per iteration so we
// never time against freed/mutated state. `make` produces fresh inputs, `run` is the timed call.
const stages = [
    {name: 'build', make: () => null, run: () => build(lng, lat, r)},
    {name: 'scan', make: () => build(lng, lat, r), run: s => scan(s)},
    {name: 'arcs', make: () => { const s = build(lng, lat, r); return [s, scan(s)]; }, run: ([s, sc]) => arcs(s, sc)},
    {name: 'stitch', make: () => { const s = build(lng, lat, r); const sc = scan(s); return [s, sc, arcs(s, sc)]; }, run: ([s, sc, a]) => stitch(s, sc, a)},
    {name: 'assemble', make: () => { const s = build(lng, lat, r); const sc = scan(s); const a = arcs(s, sc); return [s, sc, a, stitch(s, sc, a)]; }, run: ([s, sc, a, ri]) => assemble(s, sc, a, ri)},
    {name: 'sample', make: () => { const s = build(lng, lat, r); const sc = scan(s); const a = arcs(s, sc); return assemble(s, sc, a, stitch(s, sc, a)); }, run: t => sample(t)},
];

function measure({make, run}) {
    const times = [];
    for (let i = 0; i < WARMUP + RUNS; i++) {
        const input = make();
        const t0 = performance.now();
        run(input);
        const dt = performance.now() - t0;
        if (i >= WARMUP) times.push(dt);
    }
    times.sort((a, b) => a - b);
    return {min: times[0], median: times[times.length >> 1]};
}

console.log(`circle-union benchmark — ${n.toLocaleString()} circles, ${WARMUP} warmup + ${RUNS} runs\n`);
console.log(`${'stage'.padEnd(12)} ${'min'.padStart(9)} ${'median'.padStart(9)}`);

let totalMin = 0, totalMedian = 0;
for (const stage of stages) {
    const {min, median} = measure(stage);
    totalMin += min; totalMedian += median;
    console.log(`${stage.name.padEnd(12)} ${min.toFixed(2).padStart(7)}ms ${median.toFixed(2).padStart(7)}ms`);
}
console.log(`${'─'.repeat(32)}`);
console.log(`${'total'.padEnd(12)} ${totalMin.toFixed(2).padStart(7)}ms ${totalMedian.toFixed(2).padStart(7)}ms`);
