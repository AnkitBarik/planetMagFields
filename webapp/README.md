# planetMagFields 3D web app

A static, browser-only 3D explorer of planetary magnetic fields built with
[three.js](https://threejs.org). Everything is computed in the browser from the
Gauss coefficients shipped with `planetmagfields`, so it can be hosted on any
static server (it is deployed to GitHub Pages at
<https://ankitbarik.github.io/planetMagFields/app/> by
`.github/workflows/documentation.yml`).

## Run locally

```bash
python3 webapp/build_data.py   # writes webapp/data/models.js and coastlines_110m.js
```

Then just open `webapp/index.html` in a browser; no server is needed. Only
three.js is fetched from the CDN, so an internet connection is required. It can equally be served over
HTTP, e.g. `python3 -m http.server -d webapp 8000`.

## Layout

| File | Purpose |
| --- | --- |
| `build_data.py` | Exports every model's coefficients (and Earth coastlines) as JS data files using the package's own readers |
| `js/field.js` | Spherical-harmonic evaluation: Br on grids, full B vector, filtering, Lowes spectrum, dipole tilt, field-line tracing and seeding |
| `js/main.js` | three.js scene, UI, Earth time animation, hover readout, and the field-line Web Worker (started from a Blob built from `field.js`, so it also works on `file://`) |
| `js/colormaps.js`, `js/spectrum.js` | Colormaps and the SVG spectrum plot |

## Notes

* All local scripts are classic scripts (no ES modules, no `fetch`), so the
  page also works from `file://`.
* Earth's IGRF field can be animated from 1900 to 2030; while playing, the
  colour scale is fixed to the largest |Br| over the whole period.
* The ℓ/m filter (like `Planet.plot_filt`) applies to the surface, field lines,
  hover readout and spectrum. Dipole tilt and the energy-partition metrics are
  always for the unfiltered model.
* Field lines are traced with RK4 directly on the analytic potential field, so
  no 3D grid extrapolation is needed. Seeds come from a Fibonacci lattice by
  default; the |Br|-weighted strategy of the Streamlit app is also available.
* Coefficients are evaluated with Schmidt semi-normalised Legendre functions
  without the Condon-Shortley phase, at true geographic longitude.
* The view is encoded in the URL (`?planet=jupiter&model=jrm09&r=1.2`), so links
  can be shared.
