// Classic script: exposes self.PMF.drawTrajectoryPlot.
// Field components along a trajectory as one SVG line chart with a shared y
// axis, a legend and a crosshair tooltip.
(function (root) {
  'use strict';

  const NS = 'http://www.w3.org/2000/svg';
  const el = (name, attrs = {}, text) => {
    const e = document.createElementNS(NS, name);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    return e;
  };

  function fmt(v) {
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e5 || a < 1e-2)) return v.toExponential(2);
    return v.toFixed(a >= 100 ? 0 : a >= 10 ? 1 : 2);
  }

  // Tick labels: integers stay integers, no trailing zeros.
  const fmtTick = (v) => (Math.abs(v) >= 1e5 || (v !== 0 && Math.abs(v) < 1e-3)
    ? v.toExponential(1) : String(+v.toPrecision(6)));

  // "Nice" axis ticks covering [lo, hi].
  function ticks(lo, hi, n = 5) {
    if (!(hi > lo)) { hi = lo + 1; }
    const raw = (hi - lo) / n;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((f) => f * mag).find((s) => s >= raw);
    const t0 = Math.floor(lo / step) * step, t1 = Math.ceil(hi / step) * step;
    const out = [];
    for (let t = t0; t <= t1 + step * 1e-9; t += step) out.push(+t.toPrecision(12));
    return out;
  }

  // Keep at most ~2*buckets points, preserving each bucket's min and max so
  // narrow peaks survive.
  function decimate(x, y, buckets = 700) {
    const n = x.length;
    if (n <= 2 * buckets) return [...Array(n).keys()];
    const idx = [];
    const size = n / buckets;
    for (let b = 0; b < buckets; b++) {
      const i0 = Math.floor(b * size), i1 = Math.min(n, Math.floor((b + 1) * size));
      let imin = i0, imax = i0;
      for (let i = i0; i < i1; i++) {
        if (y[i] < y[imin]) imin = i;
        if (y[i] > y[imax]) imax = i;
      }
      if (imin < imax) idx.push(imin, imax); else if (imax < imin) idx.push(imax, imin); else idx.push(imin);
    }
    return idx;
  }

  /**
   * opts: { x: Float64Array, xLabel, unit,
   *         series: [{label, color, values, width, measured?}],
   *         segments?: [{start, end, label}]  (index ranges, end exclusive),
   *         onHover: (index | null) => void }
   * Model values are drawn as lines; optional measured values as dots.
   */
  function drawTrajectoryPlot(svg, legend, tip, opts) {
    const { x, xLabel, unit, series, onHover } = opts;
    // Draw at the element's real pixel width so text stays at its CSS size.
    const W = Math.round(svg.getBoundingClientRect().width) || 420;
    const H = Math.round(Math.min(420, Math.max(240, W * 0.62)));
    const mr = 14, mb = 40;
    let mt = 10;
    let ml = 58; // widened below to fit the longest y tick label
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.replaceChildren();
    const n = x.length;
    if (!n) return;

    // Legend: series names, so identity is never by colour alone.
    const hasMeasured = series.some((s) => s.measured);
    legend.replaceChildren(...series.map((s) => {
      const li = document.createElement('li');
      const sw = document.createElement('span');
      sw.className = 'sw';
      sw.style.background = s.color;
      li.append(sw, s.label);
      return li;
    }));
    if (hasMeasured) {
      const key = (cls, text) => {
        const li = document.createElement('li');
        li.className = 'mk';
        const m = document.createElement('span');
        m.className = cls;
        li.append(m, text);
        return li;
      };
      legend.append(key('mk-line', 'model'), key('mk-dot', 'measured'));
    }

    // Segments (e.g. separate perijove passes) are laid out side by side, each
    // with its own x range and ticks, so long gaps between passes take no space
    // and lines never join across them.
    const segs = (opts.segments && opts.segments.length ? opts.segments : [{ start: 0, end: n }])
      .map((sg) => {
        let a = Infinity, b = -Infinity;
        for (let i = sg.start; i < sg.end; i++) { if (x[i] < a) a = x[i]; if (x[i] > b) b = x[i]; }
        if (!(b > a)) b = a + 1;
        return { ...sg, x0: a, x1: b };
      });
    const multi = segs.length > 1;
    if (multi) mt = 26; // room for the pass labels

    let y0 = Infinity, y1 = -Infinity;
    for (const s of series) {
      for (const v of s.values) { if (v < y0) y0 = v; if (v > y1) y1 = v; }
      if (s.measured) for (const v of s.measured) { if (v < y0) y0 = v; if (v > y1) y1 = v; }
    }
    const yt = ticks(Math.min(y0, 0), Math.max(y1, 0));
    y0 = yt[0]; y1 = yt[yt.length - 1];
    ml = Math.max(ml, 30 + 7 * Math.max(...yt.map((t) => fmtTick(t).length)));
    const Y = (v) => mt + ((y1 - v) / (y1 - y0)) * (H - mt - mb);

    // Horizontal layout: width shared by duration, with a gap between segments.
    const gap = multi ? 18 : 0;
    const avail = W - ml - mr - gap * (segs.length - 1);
    const total = segs.reduce((a, sg) => a + (sg.x1 - sg.x0), 0);
    let cursor = ml;
    for (const sg of segs) {
      sg.px0 = cursor;
      sg.pw = avail * (sg.x1 - sg.x0) / total;
      cursor += sg.pw + gap;
    }
    const px = new Float64Array(n);
    for (const sg of segs) {
      for (let i = sg.start; i < sg.end; i++) px[i] = sg.px0 + ((x[i] - sg.x0) / (sg.x1 - sg.x0)) * sg.pw;
    }

    // Recessive grid and axes.
    for (const t of yt) {
      svg.append(el('line', { x1: ml, x2: W - mr, y1: Y(t), y2: Y(t), class: t === 0 ? 'zero' : 'grid' }));
      svg.append(el('text', { x: ml - 8, y: Y(t) + 4, 'text-anchor': 'end', class: 'tick' }, fmtTick(t)));
    }
    segs.forEach((sg, k) => {
      const nt = Math.max(3, Math.round((5 * sg.pw) / (W - ml - mr)));
      for (const t of ticks(sg.x0, sg.x1, nt).filter((t) => t >= sg.x0 - 1e-9 && t <= sg.x1 + 1e-9)) {
        const tx = sg.px0 + ((t - sg.x0) / (sg.x1 - sg.x0)) * sg.pw;
        svg.append(el('text', { x: tx, y: H - mb + 17, 'text-anchor': 'middle', class: 'tick' }, fmtTick(t)));
      }
      if (multi) {
        svg.append(el('text', { x: sg.px0 + sg.pw / 2, y: 14, 'text-anchor': 'middle', class: 'seglabel' },
          sg.label || `Pass ${k + 1}`));
        if (k > 0) {
          const gx = sg.px0 - gap / 2;
          svg.append(el('line', { x1: gx, x2: gx, y1: mt, y2: H - mb, class: 'segsep' }));
        }
      }
    });
    svg.append(el('text', { x: (ml + W - mr) / 2, y: H - 6, 'text-anchor': 'middle', class: 'axis' }, xLabel));
    svg.append(el('text', {
      x: 14, y: (mt + H - mb) / 2, 'text-anchor': 'middle', class: 'axis',
      transform: `rotate(-90 14 ${(mt + H - mb) / 2})`,
    }, `B (${unit})`));

    // Measured values as dots, under the model lines.
    const stride = Math.max(1, Math.ceil(n / 2500));
    for (const s of series) {
      if (!s.measured) continue;
      const g = el('g', { fill: s.color, 'fill-opacity': 0.75 });
      for (let i = 0; i < n; i += stride) {
        g.append(el('circle', { cx: px[i].toFixed(1), cy: Y(s.measured[i]).toFixed(1), r: 1.7, class: 'mdot' }));
      }
      svg.append(g);
    }

    // Model lines, one polyline per segment.
    for (const s of series) {
      for (const sg of segs) {
        const xs = x.subarray(sg.start, sg.end), vs = s.values.subarray(sg.start, sg.end);
        const idx = decimate(xs, vs, Math.max(100, Math.round(700 * sg.pw / (W - ml - mr))));
        svg.append(el('polyline', {
          class: 'line', stroke: s.color, 'stroke-width': s.width || 2,
          points: idx.map((i) => `${px[sg.start + i].toFixed(1)},${Y(vs[i]).toFixed(1)}`).join(' '),
        }));
      }
    }

    // Hover layer: crosshair, dots and a tooltip with every series' value.
    const hover = el('g', { visibility: 'hidden' });
    const cross = el('line', { y1: mt, y2: H - mb, class: 'cross' });
    hover.append(cross);
    const dots = series.map((s) => {
      const c = el('circle', { r: 4, fill: s.color, class: 'dot' });
      hover.append(c);
      return c;
    });
    svg.append(hover);
    const hit = el('rect', { x: ml, y: mt, width: W - ml - mr, height: H - mt - mb, fill: 'transparent' });
    svg.append(hit);

    // Pixel positions increase monotonically when each segment's x does;
    // otherwise fall back to a linear scan.
    let monotonic = true;
    for (let i = 1; i < n && monotonic; i++) if (px[i] < px[i - 1]) monotonic = false;
    const nearest = (pv) => {
      if (monotonic) {
        let lo = 0, hi = n - 1;
        while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (px[mid] < pv) lo = mid; else hi = mid; }
        return Math.abs(px[lo] - pv) <= Math.abs(px[hi] - pv) ? lo : hi;
      }
      let best = 0;
      for (let i = 1; i < n; i++) if (Math.abs(px[i] - pv) < Math.abs(px[best] - pv)) best = i;
      return best;
    };
    const segOf = (i) => segs.findIndex((sg) => i >= sg.start && i < sg.end);
    const wrap = svg.parentElement;
    const show = (evt) => {
      const pt = svg.createSVGPoint();
      pt.x = evt.clientX; pt.y = evt.clientY;
      const p = pt.matrixTransform(svg.getScreenCTM().inverse());
      const i = nearest(p.x);
      const cx = px[i];
      cross.setAttribute('x1', cx); cross.setAttribute('x2', cx);
      series.forEach((s, k) => { dots[k].setAttribute('cx', cx); dots[k].setAttribute('cy', Y(s.values[i])); });
      hover.setAttribute('visibility', 'visible');
      const sg = segs[segOf(i)];
      tip.innerHTML = (multi ? `<div><b>${sg.label || `Pass ${segOf(i) + 1}`}</b></div>` : '') +
        `<div>${xLabel}: <b>${fmtTick(x[i])}</b></div>` +
        (hasMeasured ? '<div style="color:var(--muted)">model / measured</div>' : '') +
        series.map((s) => `<div><span class="k" style="background:${s.color}"></span>${s.label}: <b>${fmt(s.values[i])}</b>` +
          (s.measured ? ` / <b>${fmt(s.measured[i])}</b>` : '') + ` ${unit}</div>`).join('');
      tip.hidden = false;
      const wr = wrap.getBoundingClientRect();
      const left = evt.clientX - wr.left + 14, top = evt.clientY - wr.top + 14;
      tip.style.left = `${Math.min(left, wr.width - tip.offsetWidth - 4)}px`;
      tip.style.top = `${Math.min(top, wr.height - tip.offsetHeight - 4)}px`;
      if (onHover) onHover(i);
    };
    const hide = () => {
      hover.setAttribute('visibility', 'hidden');
      tip.hidden = true;
      if (onHover) onHover(null);
    };
    hit.addEventListener('pointermove', show);
    hit.addEventListener('pointerdown', show);
    hit.addEventListener('pointerleave', hide);
  }

  root.PMF = Object.assign(root.PMF || {}, { drawTrajectoryPlot });
})(self);
