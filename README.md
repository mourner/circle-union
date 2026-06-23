# circle-union

A really fast library for computing the **union of geographic circles** as coverage polygons in JavaScript.

Given `N` points each with a radius, it computes their combined coverage area as a GeoJSON `MultiPolygon` — shells with holes — in milliseconds, even for tens of thousands of circles. Useful for things like cell-tower coverage, service-area maps, sensor ranges, or any "what's within `r` km of these points" question.

Rather than buffering each circle into a many-sided polygon and running a general-purpose boolean union, it works directly with the **arcs** that bound a union of disks — whose complexity is only `O(n)` — and does all geometry on the unit sphere. That makes it both fast and geodesic-exact: no projection distortion, and no special-casing around the antimeridian or the poles.

## Usage

```js
import {CircleUnion} from 'circle-union';

// reserve space for a known number of circles
const u = new CircleUnion(circles.length);

// add each circle (lng, lat in degrees, radius in km)
for (const {lng, lat, r} of circles) u.add(lng, lat, r);

// compute the union as a GeoJSON MultiPolygon
const geojson = u.geojson(); // {type: 'MultiPolygon', coordinates: [...]}
```

## Install

Install with NPM: `npm install circle-union`, then import as a module:

```js
import {CircleUnion} from 'circle-union';
```

Or use it directly in the browser with [jsDelivr](https://www.jsdelivr.com/esm):

```html
<script type="module">
    import {CircleUnion} from 'https://cdn.jsdelivr.net/npm/circle-union/+esm';
</script>
```

## API

#### `new CircleUnion(numItems)`

Creates a builder that will hold a given number of circles (`numItems`).

#### `u.add(lng, lat, r)`

Adds a circle centered at `lng`, `lat` (degrees) with radius `r` (km). Returns a zero-based index. Throws if you add more circles than reserved.

#### `u.geojson([options])`

Computes the union and returns it as a GeoJSON `MultiPolygon` (`{type, coordinates}`). The boundary arcs are sampled into vertices adaptively. Accepts an optional options object:

- `tolerance`: maximum arc-to-chord deviation in km (`0.005`, ≈5 m, by default) — smaller values produce smoother, denser output.
- `minPoints`: floor on the number of vertices per full circle (`24` by default), so even tiny circles stay round.

```js
const geojson = u.geojson({tolerance: 0.01, minPoints: 32});
```

#### `u.arcs()`

Returns the exact, resolution-independent arc topology behind the union, as nested GeoJSON-like arrays:

```
arc = [lng, lat, radius, startAngle, endAngle]
ring = [arc, ...]
polygon = [ring, ...]
result = [polygon, ...]
```

Use this if you want to render or measure the boundary without sampling it into line segments. Both `arcs()` and `geojson()` cache their work, so calling them repeatedly (or together) is cheap.

## Performance

Union of an [OpenCelliD](https://opencellid.org) export of cell towers over Ukraine into a GeoJSON `MultiPolygon`, on a MacBook Pro (M1 Pro, Node v24) — against general polygon-union approaches over matching 24-segment circles: [martinez-polygon-clipping](https://github.com/w8r/martinez) (a fast clipping library) and the usual `turf.circle` + `turf.union`:

circles | circle-union | martinez | turf
--: | --: | --: | --:
1,000 | 1.8 ms | 130 ms | 870 ms
4,000 | 4.8 ms | 560 ms | 9.7 s
8,000 | 6.9 ms | 930 ms | 39 s
23,467 | **18.5 ms** | 1.9 s | out of memory

Working directly with boundary arcs instead of densified polygons keeps the cost roughly linear: circle-union stays ~100× ahead of even a fast general clipper, while `turf` blows up and runs out of memory before it can finish the full set.

## Development

```sh
npm test          # invariants + golden area + independent oracle checks (~2 s, deterministic)
npm run bench     # full-pipeline timing on the real fixture
npm run preview   # static-serve the repo, then open http://localhost:8000/preview/
```

Correctness is verified against an independent brute-force membership oracle — a point is in the union iff it lies inside some disk — that shares no code with the arc algorithm, and the union area is pinned to a golden snapshot.
