// Classic script: exposes self.PMF.drawSpectrum.
(function (root) {
  'use strict';

  // Lowes spectrum as a small SVG plot with a logarithmic y axis.
  const NS = 'http://www.w3.org/2000/svg';

  function el(name, attrs = {}, text) {
    const e = document.createElementNS(NS, name);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    return e;
  }

  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const pow10 = (k) => '10' + String(k).split('').map((c) => SUP[c]).join('');

  /**
   * E: spectrum to plot (possibly filtered). Efull: optional unfiltered
   * spectrum, drawn dashed in grey behind E for comparison.
   */
  function drawSpectrum(svg, E, lmax, Efull) {
    const W = 480, H = 360, ml = 64, mr = 16, mt = 16, mb = 48;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.replaceChildren();

    const pts = (S) => {
      const ls = [], vs = [];
      if (S) for (let l = 1; l <= lmax; l++) if (S[l] > 0) { ls.push(l); vs.push(Math.log10(S[l])); }
      return { ls, vs };
    };
    const { ls, vs } = pts(E);
    const full = pts(Efull);
    const allv = vs.concat(full.vs);
    if (!allv.length) return;
    let y0 = Math.floor(Math.min(...allv)), y1 = Math.ceil(Math.max(...allv));
    if (y1 === y0) y1 = y0 + 1;
    const x = (l) => ml + ((l - 0.5) / (lmax + 0.5)) * (W - ml - mr);
    const y = (v) => mt + ((y1 - v) / (y1 - y0)) * (H - mt - mb);

    const step = Math.max(1, Math.ceil((y1 - y0) / 7));
    for (let k = y0; k <= y1; k += step) {
      svg.append(el('line', { x1: ml, x2: W - mr, y1: y(k), y2: y(k), class: 'grid' }));
      svg.append(el('text', { x: ml - 8, y: y(k) + 4, 'text-anchor': 'end', class: 'tick' }, pow10(k)));
    }
    const dl = lmax <= 15 ? 1 : 2;
    for (let l = 1; l <= lmax; l += dl) {
      svg.append(el('line', { x1: x(l), x2: x(l), y1: mt, y2: H - mb, class: 'grid' }));
      svg.append(el('text', { x: x(l), y: H - mb + 18, 'text-anchor': 'middle', class: 'tick' }, l));
    }
    svg.append(el('text', { x: (ml + W - mr) / 2, y: H - 8, 'text-anchor': 'middle', class: 'axis' }, 'ℓ'));
    svg.append(el('text', {
      x: 16, y: (mt + H - mb) / 2, 'text-anchor': 'middle', class: 'axis',
      transform: `rotate(-90 16 ${(mt + H - mb) / 2})`,
    }, 'Rℓ (nT²)'));

    if (full.ls.length) {
      svg.append(el('polyline', {
        points: full.ls.map((l, i) => `${x(l)},${y(full.vs[i])}`).join(' '), class: 'series-full',
      }));
      full.ls.forEach((l, i) => svg.append(el('circle', { cx: x(l), cy: y(full.vs[i]), r: 2.5, class: 'dot-full' })));
      svg.append(el('line', { x1: W - mr - 130, x2: W - mr - 110, y1: mt + 10, y2: mt + 10, class: 'series-full' }));
      svg.append(el('text', { x: W - mr - 104, y: mt + 14, class: 'legend' }, 'full'));
      svg.append(el('line', { x1: W - mr - 70, x2: W - mr - 50, y1: mt + 10, y2: mt + 10, class: 'series' }));
      svg.append(el('text', { x: W - mr - 44, y: mt + 14, class: 'legend' }, 'filtered'));
    }
    svg.append(el('polyline', {
      points: ls.map((l, i) => `${x(l)},${y(vs[i])}`).join(' '), class: 'series',
    }));
    ls.forEach((l, i) => {
      const c = el('circle', { cx: x(l), cy: y(vs[i]), r: 4.5, class: 'dot' });
      c.append(el('title', {}, `ℓ = ${l}: ${E[l].toExponential(3)} nT²`));
      svg.append(c);
    });
  }

  root.PMF = Object.assign(root.PMF || {}, { drawSpectrum });
})(self);
