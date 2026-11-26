#!/usr/bin/env python3
"""
Export Gauss coefficients and auxiliary data for the static three.js web app.

The web app evaluates the field entirely in the browser, so all it needs is the
coefficients of every model (plus Earth coastlines). They are written as small
JavaScript files that set globals (``self.PMF_DATA``, ``self.PMF_COAST``) rather
than JSON, so the app also works when index.html is opened from disk, where
browsers block fetch(). This script reuses the
planetmagfields readers so the app always sees exactly the same coefficients as
the Python package.

Usage
-----
    python3 webapp/build_data.py [--out webapp/data]

Coefficient convention in the exported data
-------------------------------------------
``g`` and ``h`` are flat triangular arrays indexed as ``l*(l+1)/2 + m`` (the
same layout as ``planetmagfields.libgauss.gen_idx``), in nT, for Schmidt
semi-normalised associated Legendre functions *without* the Condon-Shortley
phase, with longitude measured eastward from the prime meridian.
"""

import argparse
import json
import os

import numpy as np

from planetmagfields import get_models, planetlist, default_model
from planetmagfields.libdata import get_data
from planetmagfields.models import model_info, model_format
from planetmagfields.utils import stdDatDir

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

IGRF_YEARS = list(range(1900, 2031, 5))

# Real spacecraft trajectories offered as examples in the app. Each CSV in
# webapp/examples has columns t, r, theta, phi and optionally measured Br,
# Btheta, Bphi (nT), with '#' comment lines describing the source.
EXAMPLES = [
    {
        "id": "voyager2-uranus",
        "file": "voyager2_uranus_48s.csv",
        "name": "Voyager 2 at Uranus (1986)",
        "planet": "uranus",
        "angles": "deg",
        "units": "nT",
        "tLabel": "Hours from closest approach",
        "credit": "Voyager 2 MAG, 48 s averages; closest approach 24 Jan 1986, "
                  "17:58 UTC at 4.19 R_U. Data: Connerney, NASA PDS",
        "url": "https://doi.org/10.17189/1520033",
    },
    {
        "id": "juno-pj01-pj03",
        "file": "juno_pj01_pj03_60s.csv",
        "name": "Juno perijoves 1 and 3 (2016)",
        "planet": "jupiter",
        "angles": "deg",
        "units": "muT",
        "tLabel": "Hours from perijove",
        "passes": ["PJ1 · 27 Aug 2016", "PJ3 · 11 Dec 2016"],
        "misfitRmax": 6,
        "credit": "Juno MAG, 60 s averages, within 10 R_J of Jupiter; the first two perijoves with "
                  "science data (none was taken at PJ2). Data: Connerney, NASA PDS",
        "url": "https://doi.org/10.17189/1519711",
    },
    {
        "id": "juno-pj32-pj33",
        "file": "juno_pj32_pj33_60s.csv",
        "name": "Juno perijoves 32 and 33 (2021)",
        "planet": "jupiter",
        "angles": "deg",
        "units": "muT",
        "tLabel": "Hours from perijove",
        "passes": ["PJ32 · 21 Feb 2021", "PJ33 · 15 Apr 2021"],
        "misfitRmax": 6,
        "credit": "Juno MAG, 60 s averages, within 10 R_J of Jupiter; the last two orbits used for "
                  "JRM33. Data: Connerney, NASA PDS",
        "url": "https://doi.org/10.17189/1519711",
    },
    {
        "id": "cassini-gf-rev271-272",
        "file": "cassini_gf_rev271_272_1m.csv",
        "name": "Cassini Grand Finale, two orbits (2017)",
        "planet": "saturn",
        "angles": "deg",
        "units": "nT",
        "tLabel": "Days from first periapsis",
        "misfitRmax": 5,
        "credit": "Cassini MAG, 1 min averages; Grand Finale Revs 271 and 272 (apoapsis to apoapsis, "
                  "23 Apr – 6 May 2017), starting with the first dive between Saturn and its rings. "
                  "Shown in local time rather than longitude (the field is axisymmetric). "
                  "Data: Dougherty et al., NASA PDS",
        "url": "https://doi.org/10.17189/5rhj-sm88",
    },
]


def _round(a, nd=4):
    return [round(float(x), nd) for x in a]


def export_model(planet, model):
    info = model_info(model, planet)
    entry = {
        "id": model,
        "name": info.get("name", model),
        "url": info.get("url"),
    }

    if model_format(planet, model) == "igrf":
        # Time-dependent: store the coefficients at every 5-year epoch; the app
        # interpolates linearly in time (extrapolating outside the table, like
        # scipy's interp1d(fill_value='extrapolate') used by the package).
        gs, hs = [], []
        for year in IGRF_YEARS:
            g, h, lmax, _, mmax = get_data(stdDatDir, planet, model, year=year)
            gs.append(_round(g))
            hs.append(_round(h))
        entry.update(lmax=int(lmax), mmax=int(mmax),
                     years=IGRF_YEARS, g=gs, h=hs)
    else:
        g, h, lmax, _, mmax = get_data(stdDatDir, planet, model)
        entry.update(lmax=int(lmax), mmax=int(mmax), g=_round(g), h=_round(h))

    return entry


def write_js(path, name, obj):
    with open(path, "w") as f:
        f.write("// Generated by webapp/build_data.py - do not edit.\n")
        f.write("self.%s = " % name)
        json.dump(obj, f, separators=(",", ":"))
        f.write(";\n")


def export_coastlines(out_dir):
    src = os.path.join(HERE, "assets", "coastlines_110m.npz")
    if not os.path.exists(src):
        print("Coastline file not found, skipping:", src)
        return
    data = np.load(src)
    # Each segment is a list of [lon, lat] pairs in degrees.
    segments = [np.round(data[k], 2).tolist() for k in data.files]
    write_js(os.path.join(out_dir, "coastlines_110m.js"), "PMF_COAST", segments)


def export_examples(out_dir):
    examples = []
    for ex in EXAMPLES:
        with open(os.path.join(HERE, "examples", ex["file"])) as f:
            entry = {k: v for k, v in ex.items() if k != "file"}
            entry["csv"] = f.read()
            examples.append(entry)
    write_js(os.path.join(out_dir, "examples.js"), "PMF_EXAMPLES", examples)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=os.path.join(HERE, "data"),
                        help="Output directory (default: webapp/data)")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    planets = {}
    for planet in planetlist:
        models = sorted(str(m) for m in get_models(planet))
        planets[planet] = {
            "default_model": default_model(planet),
            "models": {m: export_model(planet, m) for m in models},
        }

    write_js(os.path.join(args.out, "models.js"), "PMF_DATA", {"planets": planets})

    export_coastlines(args.out)
    export_examples(args.out)
    print("Wrote web app data to", args.out)


if __name__ == "__main__":
    main()
