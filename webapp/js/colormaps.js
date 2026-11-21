// Classic script: exposes self.PMF.{COLORMAPS, PLASMA, sample, lut}.
(function (root) {
  'use strict';

  // Diverging colormaps, listed from the low (negative) end to the high end.
  // Control points from ColorBrewer / matplotlib.
  const hex = (s) => [1, 3, 5].map((i) => parseInt(s.slice(i, i + 2), 16));

  const brewer = (list) => list.split(' ').map(hex);

  const COLORMAPS = {
    // Blue = negative, red = positive (the Streamlit app's default look).
    RdBu: brewer('#053061 #2166ac #4393c3 #92c5de #d1e5f0 #f7f7f7 #fddbc7 #f4a582 #d6604d #b2182b #67001f'),
    seismic: [[0, 0, 77], [0, 0, 255], [255, 255, 255], [255, 0, 0], [128, 0, 0]],
    coolwarm: brewer('#3b4cc0 #6688ee #88bbff #b8d0f9 #dddddd #f5c4ad #f49a7b #e36a53 #b40426'),
    RdYlBu: brewer('#313695 #4575b4 #74add1 #abd9e9 #e0f3f8 #ffffbf #fee090 #fdae61 #f46d43 #d73027 #a50026'),
    PuOr: brewer('#2d004b #542788 #8073ac #b2abd2 #d8daeb #f7f7f7 #fee0b6 #fdb863 #e08214 #b35806 #7f3b08'),
    BrBG: brewer('#543005 #8c510a #bf812d #dfc27d #f6e8c3 #f5f5f5 #c7eae5 #80cdc1 #35978f #01665e #003c30'),
    PiYG: brewer('#276419 #4d9221 #7fbc41 #b8e186 #e6f5d0 #f7f7f7 #fde0ef #f1b6da #de77ae #c51b7d #8e0152'),
    PRGn: brewer('#00441b #1b7837 #5aae61 #a6dba0 #d9f0d3 #f7f7f7 #e7d4e8 #c2a5cf #9970ab #762a83 #40004b'),
    Spectral: brewer('#5e4fa2 #3288bd #66c2a5 #abdda4 #e6f598 #ffffbf #fee08b #fdae61 #f46d43 #d53e4f #9e0142'),
  };

  // Sequential map used to colour field lines in "rainbow" mode.
  const PLASMA = brewer('#0d0887 #46039f #7201a8 #9c179e #bd3786 #d8576b #ed7953 #fb9f3a #fdca26 #f0f921');

  /** Colour at t in [0, 1] as [r, g, b] in 0..255. */
  function sample(stops, t) {
    const x = Math.min(1, Math.max(0, t)) * (stops.length - 1);
    const i = Math.min(stops.length - 2, Math.floor(x));
    const f = x - i;
    const a = stops[i], b = stops[i + 1];
    return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
  }

  /** 256-entry lookup table (Uint8Array RGB). */
  function lut(name, invert = false) {
    const stops = COLORMAPS[name] || COLORMAPS.RdBu;
    const out = new Uint8Array(256 * 3);
    for (let k = 0; k < 256; k++) {
      const c = sample(stops, invert ? 1 - k / 255 : k / 255);
      out.set(c.map(Math.round), 3 * k);
    }
    return out;
  }

  root.PMF = Object.assign(root.PMF || {}, { COLORMAPS, PLASMA, sample, lut });
})(self);
