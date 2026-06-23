// Live preview — parses the raw OpenCelliD tower CSV in the browser, draws each tower as a circle, and
// computes the union boundary end to end with the actual library. Loaded as a module by index.html;
// `maplibregl` comes from the CDN <script> there, bare imports inside `index.js` resolve via the import map
// there. Serve the repo root and open /preview/.
/* global maplibregl */
import {CircleUnion} from '../index.js';

const R = 6371; // km
const CELLS = '../test/fixtures/ukraine-cell-id.csv';

// build a geodesic circle ring (lng/lat) with `steps` segments
function circleRing(lng, lat, radiusKm, steps) {
    const ring = [];
    const lat1 = lat * Math.PI / 180, lng1 = lng * Math.PI / 180;
    const d = radiusKm / R;
    const sinLat1 = Math.sin(lat1), cosLat1 = Math.cos(lat1);
    const sinD = Math.sin(d), cosD = Math.cos(d);
    for (let i = 0; i <= steps; i++) {
        const brng = 2 * Math.PI * i / steps;
        const lat2 = Math.asin(sinLat1 * cosD + cosLat1 * sinD * Math.cos(brng));
        const lng2 = lng1 + Math.atan2(
            Math.sin(brng) * sinD * cosLat1,
            cosD - sinLat1 * Math.sin(lat2));
        ring.push([lng2 * 180 / Math.PI, lat2 * 180 / Math.PI]);
    }
    return ring;
}

const map = new maplibregl.Map({
    container: 'map',
    style: 'https://tiles.openfreemap.org/styles/positron',
    center: [31, 49],
    zoom: 5
});
map.addControl(new maplibregl.NavigationControl());

async function load() {
    // Parse the raw tower CSV: draw each as a circle, and collect SoA arrays for the union.
    const text = await (await fetch(CELLS)).text();
    const lines = text.split('\n');
    const features = [];
    const lng = [], lat = [], r = [];
    for (const line of lines) {
        if (!line) continue;
        const c = line.split(','); // lon,lat,range_m
        const clng = +c[0], clat = +c[1], radiusKm = +c[2] / 1000;
        lng.push(clng); lat.push(clat); r.push(radiusKm);
        // fewer segments for tiny circles, more for big ones
        const steps = radiusKm > 10 ? 48 : radiusKm > 2 ? 24 : 12;
        features.push({
            type: 'Feature',
            properties: {r: radiusKm},
            geometry: {type: 'Polygon', coordinates: [circleRing(clng, clat, radiusKm, steps)]}
        });
    }

    map.addSource('circles', {type: 'geojson', data: {type: 'FeatureCollection', features}});
    map.addLayer({
        id: 'fill',
        type: 'fill',
        source: 'circles',
        paint: {'fill-color': '#1e88e5', 'fill-opacity': 0.15}
    });
    map.addLayer({
        id: 'outline',
        type: 'line',
        source: 'circles',
        paint: {'line-color': '#1565c0', 'line-width': 0.4, 'line-opacity': 0.5}
    });

    // Compute the union boundary right here, end to end, from the same circles.
    const t0 = performance.now();
    const u = new CircleUnion(lng.length);
    for (let i = 0; i < lng.length; i++) u.add(lng[i], lat[i], r[i]);
    const geometry = u.geojson();
    const ms = performance.now() - t0;

    let vertices = 0;
    for (const poly of geometry.coordinates) for (const ring of poly) vertices += ring.length;

    map.addSource('union', {type: 'geojson', data: {type: 'Feature', properties: {}, geometry}});
    map.addLayer({
        id: 'union-line',
        type: 'line',
        source: 'union',
        paint: {'line-color': '#e53935', 'line-width': 1.5}
    });

    document.getElementById('info').textContent =
        `${features.length.toLocaleString()} towers → ${geometry.coordinates.length} polygons, ` +
        `${vertices.toLocaleString()} vertices in ${ms.toFixed(0)} ms (red = union boundary)`;

    const opacity = document.getElementById('opacity');
    opacity.oninput = () => map.setPaintProperty('fill', 'fill-opacity', +opacity.value);

    const stroke = document.getElementById('stroke');
    stroke.onchange = () => map.setLayoutProperty('outline', 'visibility', stroke.checked ? 'visible' : 'none');
}

map.on('load', () => load().catch((e) => {
    document.getElementById('info').textContent = `error: ${e.message}`;
}));
