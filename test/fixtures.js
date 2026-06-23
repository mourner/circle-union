// Shared sample data for the test suite and benchmark (dev only — not shipped).
import {readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';

// The OpenCelliD Ukraine cell-tower sample: ~23k geographic disks, the workload the
// library was built for. CSV columns: radio,mcc,net,area,cell,unit,lon,lat,range,…
export const CELL_ID = fileURLToPath(new URL('./fixtures/ukraine-cell-id.csv', import.meta.url));

/** Parse a cell-tower CSV into SoA typed arrays {n, lng, lat, r} (r in km). */
export function loadCells(path = CELL_ID) {
    const lines = readFileSync(path, 'utf8').split('\n');
    const n = lines.length;
    const lng = new Float64Array(n), lat = new Float64Array(n), r = new Float64Array(n);
    let count = 0;
    for (let i = 0; i < n; i++) {
        const line = lines[i];
        if (!line) continue;
        const c = line.split(',');
        const rangeM = +c[8];
        if (!rangeM) continue;
        lng[count] = +c[6];
        lat[count] = +c[7];
        r[count] = rangeM / 1000; // m → km
        count++;
    }
    return {n: count, lng: lng.subarray(0, count), lat: lat.subarray(0, count), r: r.subarray(0, count)};
}
