// Classic script (not an ES module) so the app works both when hosted and when
// index.html is opened straight from disk (file://), where browsers block local
// module files, fetch() and workers. three.js is pulled from the CDN with
// dynamic import(), which the import map in index.html resolves.
(async () => {
const THREE = await import('three');
const [{ OrbitControls }, { Line2 }, { LineGeometry }, { LineMaterial }] = await Promise.all([
  import('three/addons/controls/OrbitControls.js'),
  import('three/addons/lines/Line2.js'),
  import('three/addons/lines/LineGeometry.js'),
  import('three/addons/lines/LineMaterial.js'),
]);
const { FieldModel, seedPoints, fibonacciSeeds, traceLine, COLORMAPS, PLASMA, lut, sample, drawSpectrum,
  drawTrajectoryPlot } = self.PMF;

const $ = (id) => document.getElementById(id);
const UNITS = { nT: [1, 'nT'], muT: [1e-3, 'μT'], Gauss: [1e-5, 'G'] };
const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1);

// ---------------------------------------------------------------------------
// Data
// ---------------------------------------------------------------------------
// Loaded by <script> tags (data/models.js, data/coastlines_110m.js).
const DATA = self.PMF_DATA;
const COAST = self.PMF_COAST || [];
if (!DATA) {
  throw new LoadError('Model data not found',
    '<p>Could not load <code>data/models.js</code>.</p>' +
    '<p>Generate it from the repository root with <code>python3 webapp/build_data.py</code>, then reload.</p>');
}

// ---------------------------------------------------------------------------
// three.js scene (z up, matching planetmagfields' spherical coordinates)
// ---------------------------------------------------------------------------
const canvas = $('scene');
const viewport = $('viewport');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x141414);

const camera = new THREE.PerspectiveCamera(35, 1, 0.01, 500);
camera.up.set(0, 0, 1);
scene.add(camera);
scene.add(new THREE.AmbientLight(0xffffff, 0.8));
const headlight = new THREE.DirectionalLight(0xffffff, 0.45);
headlight.position.set(0.4, 0.6, 1);
camera.add(headlight);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.12;
controls.addEventListener('change', requestRender);

// Unit sphere with equirectangular UVs: u = phi / 2pi, v = theta / pi.
function sphereGeometry(nphi, ntheta) {
  const pos = [], uv = [], index = [];
  for (let j = 0; j <= ntheta; j++) {
    const th = (j / ntheta) * Math.PI;
    for (let i = 0; i <= nphi; i++) {
      const ph = (i / nphi) * 2 * Math.PI;
      pos.push(Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th));
      uv.push(i / nphi, j / ntheta);
    }
  }
  for (let j = 0; j < ntheta; j++) {
    for (let i = 0; i < nphi; i++) {
      const a = j * (nphi + 1) + i, b = a + nphi + 1;
      index.push(a, b, a + 1, b, b + 1, a + 1);
    }
  }
  const geo = new THREE.BufferGeometry();
  geo.setIndex(index);
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  geo.setAttribute('normal', new THREE.Float32BufferAttribute(pos, 3));
  geo.setAttribute('uv', new THREE.Float32BufferAttribute(uv, 2));
  return geo;
}

const surfaceMat = new THREE.MeshLambertMaterial({ color: 0xffffff });
const surface = new THREE.Mesh(sphereGeometry(256, 128), surfaceMat);
scene.add(surface);
let texture = null;

const coastGroup = new THREE.Group();
scene.add(coastGroup);
const coastMat = new THREE.LineBasicMaterial({ color: 0xffffff });
{
  const pts = [];
  const toXYZ = ([lon, lat]) => {
    const th = Math.PI / 2 - (lat * Math.PI) / 180;
    const ph = (lon * Math.PI) / 180; // true geographic longitude
    return [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
  };
  for (const seg of COAST) {
    for (let k = 0; k + 1 < seg.length; k++) pts.push(...toXYZ(seg[k]), ...toXYZ(seg[k + 1]));
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
  coastGroup.add(new THREE.LineSegments(geo, coastMat));
}

const linesGroup = new THREE.Group();
scene.add(linesGroup);
const lineMat = new LineMaterial({ color: 0xffffff, linewidth: 3, vertexColors: true });

const axesGroup = new THREE.Group();
scene.add(axesGroup);
function buildAxes(limit) {
  axesGroup.clear();
  const colors = [0xff5555, 0x55dd55, 0x5599ff];
  ['x', 'y', 'z'].forEach((name, k) => {
    const dir = new THREE.Vector3().setComponent(k, 1);
    const geo = new THREE.BufferGeometry().setFromPoints([dir.clone().multiplyScalar(-limit), dir.clone().multiplyScalar(limit)]);
    axesGroup.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ color: colors[k] })));
    const c = document.createElement('canvas');
    c.width = c.height = 64;
    const ctx = c.getContext('2d');
    ctx.fillStyle = '#' + colors[k].toString(16).padStart(6, '0');
    ctx.font = 'bold 48px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(name, 32, 34);
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), depthTest: false }));
    sprite.position.copy(dir.multiplyScalar(limit * 1.06));
    sprite.scale.setScalar(limit * 0.08);
    axesGroup.add(sprite);
  });
}

// ---------------------------------------------------------------------------
// State and controls
// ---------------------------------------------------------------------------
const ui = {};
for (const id of ['planet', 'model', 'year', 'r', 'units', 'resolution', 'show-spectrum', 'show-axes',
  'show-coast', 'coast-color', 'cmap', 'invert', 'sym', 'sat', 'vmin', 'vmax', 'show-lines', 'nlines',
  'seeding', 'rmax', 'line-color', 'lw', 'lmin', 'lmaxf', 'mmin', 'mmaxf', 'play', 'speed', 'traj-angles']) ui[id] = $(id);

const val = {
  get planet() { return ui.planet.value; },
  get model() { return ui.model.value; },
  get year() { return this.planet === 'earth' ? parseFloat(ui.year.value) : null; },
  get r() { return parseFloat(ui.r.value); },
  get nphi() { return parseInt(ui.resolution.value, 10); },
  get unit() { return UNITS[ui.units.value]; },
  get rmax() { return parseFloat(ui.rmax.value); },
};

for (const p of Object.keys(DATA.planets)) ui.planet.add(new Option(cap(p), p));
for (const c of Object.keys(COLORMAPS)) ui.cmap.add(new Option(c, c));
ui.year.value = Math.min(2030, new Date().getFullYear());

function fillModels(preferred) {
  const P = DATA.planets[val.planet];
  ui.model.replaceChildren();
  for (const [id, m] of Object.entries(P.models)) ui.model.add(new Option(m.name, id));
  ui.model.value = preferred && P.models[preferred] ? preferred : P.default_model;
}

// Restore state from the URL so views can be shared.
const params = new URLSearchParams(location.search);
ui.planet.value = DATA.planets[params.get('planet')] ? params.get('planet') : 'earth';
fillModels(params.get('model'));
for (const k of ['year', 'r', 'rmax', 'nlines']) if (params.has(k)) ui[k].value = params.get(k);
if (params.has('units') && UNITS[params.get('units')]) ui.units.value = params.get('units');
if (params.has('cmap') && COLORMAPS[params.get('cmap')]) ui.cmap.value = params.get('cmap');
if (params.has('seeding')) ui.seeding.value = params.get('seeding');
if (params.has('sat')) ui.sat.value = params.get('sat');

// Colormap picker: shows a preview of every colormap (following "Invert") and
// drives the hidden <select id="cmap">, so the rest of the app reads ui.cmap.
// Arrow keys preview colormaps live on the globe; Escape restores the original.
const cmapBtn = $('cmap-btn');
const cmapList = $('cmap-list');
let cmapOriginal = null;
for (const name of Object.keys(COLORMAPS)) {
  const li = document.createElement('li');
  li.id = `cmap-opt-${name}`;
  li.setAttribute('role', 'option');
  li.dataset.name = name;
  li.tabIndex = -1;
  li.innerHTML = '<span class="cmap-swatch" aria-hidden="true"></span>' +
    `<span class="cmap-name">${name}</span><span class="check" aria-hidden="true">✓</span>`;
  cmapList.append(li);
}
function gradientCSS(name) {
  const t = lut(name, ui.invert.checked);
  const stops = [];
  for (let k = 0; k <= 12; k++) {
    const i = Math.round((k / 12) * 255);
    stops.push(`rgb(${t[3 * i]},${t[3 * i + 1]},${t[3 * i + 2]}) ${((k / 12) * 100).toFixed(1)}%`);
  }
  return `linear-gradient(to right, ${stops.join(',')})`;
}
function refreshCmapPicker() {
  const cur = ui.cmap.value;
  cmapBtn.querySelector('.cmap-swatch').style.background = gradientCSS(cur);
  cmapBtn.querySelector('.cmap-name').textContent = cur;
  for (const li of cmapList.children) {
    li.querySelector('.cmap-swatch').style.background = gradientCSS(li.dataset.name);
    li.setAttribute('aria-selected', String(li.dataset.name === cur));
  }
}
function applyCmap(name) {
  if (name !== ui.cmap.value) {
    ui.cmap.value = name;
    ui.cmap.dispatchEvent(new Event('input'));
  }
  refreshCmapPicker();
}
function openCmapList() {
  cmapOriginal = ui.cmap.value;
  cmapList.hidden = false;
  cmapBtn.setAttribute('aria-expanded', 'true');
  const sel = cmapList.querySelector('[aria-selected="true"]') || cmapList.firstElementChild;
  sel.focus({ preventScroll: true });
  sel.scrollIntoView({ block: 'nearest' });
}
function closeCmapList(focusButton) {
  if (cmapList.hidden) return;
  cmapList.hidden = true;
  cmapBtn.setAttribute('aria-expanded', 'false');
  for (const li of cmapList.children) li.classList.remove('active');
  if (focusButton) cmapBtn.focus();
}
cmapBtn.addEventListener('click', () => (cmapList.hidden ? openCmapList() : closeCmapList(false)));
cmapBtn.addEventListener('keydown', (e) => {
  if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); openCmapList(); }
});
cmapList.addEventListener('click', (e) => {
  const li = e.target.closest('li');
  if (!li) return;
  applyCmap(li.dataset.name);
  closeCmapList(true);
});
cmapList.addEventListener('focusin', (e) => {
  const li = e.target.closest('li');
  if (!li) return;
  for (const o of cmapList.children) o.classList.toggle('active', o === li);
  cmapList.setAttribute('aria-activedescendant', li.id);
});
cmapList.addEventListener('keydown', (e) => {
  const items = [...cmapList.children];
  let i = items.indexOf(document.activeElement);
  switch (e.key) {
    case 'ArrowDown': i = Math.min(items.length - 1, i + 1); break;
    case 'ArrowUp': i = Math.max(0, i - 1); break;
    case 'Home': i = 0; break;
    case 'End': i = items.length - 1; break;
    case 'Enter':
    case ' ':
      e.preventDefault();
      applyCmap(items[i].dataset.name);
      closeCmapList(true);
      return;
    case 'Escape':
      e.preventDefault();
      applyCmap(cmapOriginal);
      closeCmapList(true);
      return;
    case 'Tab':
      closeCmapList(false);
      return;
    default:
      return;
  }
  e.preventDefault();
  items[i].focus();
  items[i].scrollIntoView({ block: 'nearest' });
  applyCmap(items[i].dataset.name); // live preview on the globe
});
document.addEventListener('pointerdown', (e) => {
  if (!cmapList.hidden && !e.target.closest('.cmap-picker')) closeCmapList(false);
});
ui.invert.addEventListener('input', refreshCmapPicker);
refreshCmapPicker();

// Spherical-harmonic filter: slider bounds follow the selected model.
function resetFilter() {
  const m = DATA.planets[val.planet].models[val.model];
  for (const id of ['lmin', 'lmaxf']) { ui[id].max = m.lmax; }
  for (const id of ['mmin', 'mmaxf']) { ui[id].max = m.lmax; }
  ui.lmin.value = 1;
  ui.lmaxf.value = m.lmax;
  ui.mmin.value = 0;
  ui.mmaxf.value = m.lmax;
  $('m-filter').hidden = m.mmax === 0;
}
const filter = () => ({
  lmin: +ui.lmin.value, lmax: +ui.lmaxf.value, mmin: +ui.mmin.value, mmax: +ui.mmaxf.value,
});
function filterActive() {
  const f = filter(), m = DATA.planets[val.planet].models[val.model];
  return f.lmin > 1 || f.lmax < m.lmax || (m.mmax > 0 && (f.mmin > 0 || f.mmax < m.lmax));
}
resetFilter();
for (const [k, id] of [['lmin', 'lmin'], ['lmax', 'lmaxf'], ['mmin', 'mmin'], ['mmax', 'mmaxf']]) {
  if (params.has(k)) ui[id].value = params.get(k);
}

function syncURL() {
  const q = new URLSearchParams({ planet: val.planet, model: val.model, r: ui.r.value });
  if (val.year !== null) q.set('year', ui.year.value);
  if (ui.units.value !== 'muT') q.set('units', ui.units.value);
  if (ui.cmap.value !== 'RdBu') q.set('cmap', ui.cmap.value);
  if (traj && (traj.source === 'synthetic' || EXAMPLES.some((e) => e.id === traj.source))) q.set('traj', traj.source);
  if (ui.sym.checked && +ui.sat.value > 0) q.set('sat', ui.sat.value);
  if (ui['show-lines'].checked) {
    q.set('nlines', ui.nlines.value);
    q.set('rmax', ui.rmax.value);
    if (ui.seeding.value !== 'fibonacci') q.set('seeding', ui.seeding.value);
  }
  if (filterActive()) {
    const f = filter();
    q.set('lmin', f.lmin); q.set('lmax', f.lmax);
    if (!$('m-filter').hidden) { q.set('mmin', f.mmin); q.set('mmax', f.mmax); }
  }
  history.replaceState(null, '', '?' + q);
}

// ---------------------------------------------------------------------------
// Pipeline: model -> surface -> colours; model -> field lines; model -> spectrum
// ---------------------------------------------------------------------------
let fullModel = null; // as published
let model = null; // after the spherical-harmonic filter; drives everything drawn
let brGrid = null; // nT, ntheta x nphi
let range = [-1, 1];
let traced = []; // Float32Array xyz per line

function buildModel() {
  const m = DATA.planets[val.planet].models[val.model];
  fullModel = m.years
    ? FieldModel.fromEpochs(m.lmax, m.years, m.g, m.h, val.year, m.mmax)
    : new FieldModel(m.lmax, m.g, m.h, m.mmax);
  applyFilter();
}

function applyFilter() {
  const f = filter();
  model = filterActive() ? fullModel.filtered(f.lmin, f.lmax, f.mmin, f.mmax) : fullModel;
}

function updateSurface() {
  const nphi = val.nphi, ntheta = nphi / 2;
  brGrid = model.gridBr(val.r, nphi, ntheta);
  if (!texture || texture.image.width !== nphi) {
    texture?.dispose();
    texture = new THREE.DataTexture(new Uint8Array(nphi * ntheta * 4), nphi, ntheta, THREE.RGBAFormat);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.wrapS = THREE.RepeatWrapping;
    texture.magFilter = texture.minFilter = THREE.LinearFilter;
    surfaceMat.map = texture;
    surfaceMat.needsUpdate = true;
  }
  const r = val.r;
  surface.scale.setScalar(r);
  coastGroup.scale.setScalar(r * 1.002);
  recolor();
}

let animMaxCache = { key: '', value: 0 };
function animAbsMax() {
  const f = filter();
  const key = `${val.model}|${val.r}|${f.lmin}|${f.lmax}|${f.mmin}|${f.mmax}`;
  if (animMaxCache.key !== key) {
    const m = DATA.planets[val.planet].models[val.model];
    let mx = 0;
    for (let y = +ui.year.min; y <= +ui.year.max; y += 5) {
      let fm = FieldModel.fromEpochs(m.lmax, m.years, m.g, m.h, y, m.mmax);
      if (filterActive()) fm = fm.filtered(f.lmin, f.lmax, f.mmin, f.mmax);
      for (const v of fm.gridBr(val.r, 128, 64)) mx = Math.max(mx, Math.abs(v));
    }
    animMaxCache = { key, value: mx };
  }
  return animMaxCache.value;
}

function recolor() {
  const [fac, label] = val.unit;
  let absmax = 0, dmin = Infinity, dmax = -Infinity;
  for (const v of brGrid) {
    if (v < dmin) dmin = v;
    if (v > dmax) dmax = v;
  }
  dmin *= fac;
  dmax *= fac;
  absmax = Math.max(-dmin, dmax);
  if (ui.sym.checked) {
    // While animating, use the largest |Br| over the whole time span so the
    // colours mean the same thing in every frame.
    if (playing) absmax = animAbsMax() * fac;
    // Saturation: colours saturate at (1 - s) of the maximum.
    absmax *= 1 - parseFloat(ui.sat.value);
    range = [-absmax, absmax];
    $('sat-out').textContent = +ui.sat.value > 0 ? `${Math.round(100 * ui.sat.value)}% · ±${fmt(absmax)} ${label}` : '0%';
  } else {
    const lo = parseFloat(ui.vmin.value), hi = parseFloat(ui.vmax.value);
    range = hi > lo ? [lo, hi] : [-absmax, absmax];
  }
  const table = lut(ui.cmap.value, ui.invert.checked);
  const data = texture.image.data;
  const span = range[1] - range[0] || 1;
  for (let n = 0; n < brGrid.length; n++) {
    const k = Math.max(0, Math.min(255, Math.round(((brGrid[n] * fac - range[0]) / span) * 255)));
    data[4 * n] = table[3 * k];
    data[4 * n + 1] = table[3 * k + 1];
    data[4 * n + 2] = table[3 * k + 2];
    data[4 * n + 3] = 255;
  }
  texture.needsUpdate = true;
  // Mark the colorbar ends when the field goes beyond the colour range.
  drawColorbar(table, label, dmin < range[0] - 1e-9 * span, dmax > range[1] + 1e-9 * span);
  colorTrajectory(); // follows colormap / invert changes
  requestRender();
}

function fmt(v) {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e4 || a < 1e-2)) return v.toExponential(2);
  return v.toFixed(a >= 100 ? 0 : a >= 10 ? 1 : 2);
}

function drawColorbar(table, label, clipLow = false, clipHigh = false) {
  const stops = [];
  for (let k = 0; k <= 16; k++) {
    const i = Math.round((k / 16) * 255);
    stops.push(`rgb(${table[3 * i]},${table[3 * i + 1]},${table[3 * i + 2]}) ${(k / 16) * 100}%`);
  }
  $('cb-grad').style.background = `linear-gradient(to top, ${stops.join(',')})`;
  $('cb-title').textContent = `Bᵣ (${label})`;
  const ticks = $('cb-ticks');
  ticks.replaceChildren();
  for (let k = 0; k <= 4; k++) {
    const d = document.createElement('div');
    d.style.top = `${100 - k * 25}%`;
    const prefix = k === 4 && clipHigh ? '≥ ' : k === 0 && clipLow ? '≤ ' : '';
    d.textContent = prefix + fmt(range[0] + (k / 4) * (range[1] - range[0]));
    ticks.append(d);
  }
}

// Field lines are traced in a worker built from a Blob (field.js source plus
// the handler below), which also works on file:// pages. Only the newest
// request is kept while the worker is busy, so animation and slider drags
// never build up a backlog.
function workerMain() {
  self.onmessage = (e) => {
    const { FieldModel, traceLine } = self.PMF;
    const { job, key, lmax, mmax, g, h, seeds, rSurf, rMax } = e.data;
    const model = new FieldModel(lmax, g, h, mmax);
    const lines = [];
    for (const [theta, phi] of seeds) {
      const pts = traceLine(model, theta, phi, rSurf, rMax);
      if (pts) lines.push(pts);
    }
    self.postMessage({ job, key, lines }, lines.map((l) => l.buffer));
  };
}
let worker = null;
try {
  const src = `(${self.PMF.fieldSource})(self);\n(${workerMain})();`;
  worker = new Worker(URL.createObjectURL(new Blob([src], { type: 'text/javascript' })));
} catch { /* fall back to the main thread */ }
let job = 0;
let busy = false;
let queued = null;
let shownJob = 0;
// Results for the same planet/model may be shown even if slightly stale (smooth
// animation); results for a different planet/model never are.
const sceneKey = () => `${val.planet}|${val.model}`;

function postTrace(msg) {
  if (busy) { queued = msg; return; }
  busy = true;
  worker.postMessage(msg);
}
if (worker) worker.onmessage = (e) => {
  busy = false;
  if (queued) { const m = queued; queued = null; postTrace(m); }
  const { job: j, key, lines } = e.data;
  if (j < shownJob || (j !== job && key !== sceneKey())) return;
  shownJob = j;
  traced = lines;
  if (j === job) setStatus(null);
  buildLines();
};
if (worker) worker.onerror = (e) => {
  e.preventDefault?.();
  worker = null;
  updateLines();
};

function updateLines() {
  job++;
  queued = null;
  const r = val.r, rmax = val.rmax;
  const clear = () => { traced = []; shownJob = job; buildLines(); };
  if (!ui['show-lines'].checked) { clear(); setStatus(null); return; }
  if (r >= rmax) {
    clear();
    setStatus(`Surface radius (${r.toFixed(2)} Rₚ) ≥ max field line radius (${rmax.toFixed(1)} Rₚ): field lines hidden.`);
    return;
  }
  const n = parseInt(ui.nlines.value, 10);
  let seeds;
  if (ui.seeding.value === 'weighted') {
    seeds = seedPoints(model.gridBr(r, 128, 64), 128, 64, n);
  } else {
    seeds = fibonacciSeeds(n);
  }
  if (!worker) {
    traced = seeds.map(([th, ph]) => traceLine(model, th, ph, r, rmax)).filter(Boolean);
    shownJob = job;
    setStatus(null);
    buildLines();
    return;
  }
  if (!playing) setStatus('Tracing field lines…');
  postTrace({ job, key: sceneKey(), lmax: model.lmax, mmax: model.mmax, g: model.g, h: model.h, seeds, rSurf: r, rMax: rmax });
}

function buildLines() {
  for (const obj of linesGroup.children) obj.geometry.dispose();
  linesGroup.clear();
  const choice = ui['line-color'].value;
  lineMat.linewidth = parseFloat(ui.lw.value);
  traced.forEach((pts, i) => {
    const col = choice === 'rainbow'
      ? new THREE.Color().setRGB(...sample(PLASMA, traced.length > 1 ? i / (traced.length - 1) : 0).map((c) => c / 255), THREE.SRGBColorSpace)
      : new THREE.Color(choice);
    const colors = new Float32Array(pts.length);
    for (let k = 0; k < pts.length; k += 3) colors.set([col.r, col.g, col.b], k);
    const geo = new LineGeometry();
    geo.setPositions(pts);
    geo.setColors(colors);
    const line = new Line2(geo, lineMat);
    line.computeLineDistances();
    linesGroup.add(line);
  });
  requestRender();
}

function updateSpectrum() {
  if (!ui['show-spectrum'].checked) return;
  const s = fullModel.spectrum(val.r);
  if (model !== fullModel) drawSpectrum($('spectrum'), model.spectrum(val.r).E, model.lmax, s.E);
  else drawSpectrum($('spectrum'), s.E, model.lmax);
  const pct = (x) => `${(100 * x).toFixed(1)}%`;
  $('s-dip').textContent = pct(s.dipolarity);
  $('s-diptot').textContent = pct(s.dipTot);
  $('s-axi').textContent = pct(s.axi / s.tot);
  $('s-symm').textContent = pct(s.symm / s.tot);
}

// ---------------------------------------------------------------------------
// Trajectory: the field along a user-supplied path (r, theta, phi), drawn on
// the globe and plotted in the side panel. Data never leaves the browser.
// ---------------------------------------------------------------------------
// Components use the first three dark-mode categorical slots (validated for
// colour-blind separation on this surface); |B|, a magnitude, is neutral ink.
const TRAJ_SERIES = [
  { key: 'br', label: 'Bᵣ', color: '#3987e5' },
  { key: 'bt', label: 'Bθ', color: '#d95926' },
  { key: 'bp', label: 'Bφ', color: '#199e70' },
  { key: 'b', label: '|B|', color: '#e6e6e6', width: 2.5 },
];
const TRAJ_COLOR = 0xf5b400;
let traj = null;   // { r, th, ph (as in file), t, tLabel, dropped, source, unit, theta, phi (radians), xyz, rmaxPath }
let trajB = null;  // { br, bt, bp, b } in nT
const trajGroup = new THREE.Group();
scene.add(trajGroup);
// The path is coloured by Br / max|Br| along it with the surface colormap,
// over a dark outline so it stays visible against the coloured globe.
const TRAJ_WIDTH = 3.5;
const trajMat = new LineMaterial({ color: 0xffffff, linewidth: TRAJ_WIDTH, vertexColors: true });
// No depth writes: the outline can't cover the coloured line drawn after it,
// but the globe still hides it where the path goes behind the planet.
const trajOutlineMat = new LineMaterial({ color: 0x0b0b0b, linewidth: TRAJ_WIDTH + 3, depthWrite: false });
const markerGeo = new THREE.SphereGeometry(1, 20, 14);
const trajStart = new THREE.Mesh(markerGeo, new THREE.MeshBasicMaterial({ color: TRAJ_COLOR }));
const trajHover = new THREE.Mesh(markerGeo, new THREE.MeshBasicMaterial({ color: 0xffffff }));
trajStart.visible = trajHover.visible = false;
scene.add(trajStart, trajHover);

const minmax = (a) => {
  let lo = Infinity, hi = -Infinity;
  for (const v of a) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return [lo, hi];
};

function parseTrajectoryCSV(text) {
  const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !l.startsWith('#'));
  if (!lines.length) throw new Error('No data found.');
  const split = (l) => l.split(/[\s,;]+/).filter(Boolean).map((c) => c.replace(/^["']|["']$/g, ''));
  const isNum = (c) => c !== '' && Number.isFinite(Number(c));
  let header = split(lines[0]);
  let start = 1;
  if (header.every(isNum)) { header = ['r', 'theta', 'phi', 't'].slice(0, header.length); start = 0; }
  const names = header.map((h) => h.toLowerCase());
  const col = (...alts) => names.findIndex((n) => alts.includes(n));
  const ir = col('r', 'radius'), it = col('theta', 'colat', 'colatitude'), ip = col('phi', 'lon', 'longitude');
  const itime = col('t', 'time');
  const ipass = col('pass', 'segment');
  // Optional measured field (nT): plotted as dots and used to colour the path.
  const im = [col('br', 'b_r'), col('btheta', 'bt', 'b_theta'), col('bphi', 'bp', 'b_phi')];
  const hasM = im.every((k) => k >= 0);
  if (ir < 0 || it < 0 || ip < 0) throw new Error(`Need columns r, theta and phi; found: ${header.join(', ')}.`);
  const r = [], th = [], ph = [], t = [], mbr = [], mbt = [], mbp = [], pass = [];
  let dropped = 0;
  for (let k = start; k < lines.length; k++) {
    const c = split(lines[k]);
    const v = [c[ir], c[it], c[ip]].map(Number);
    if (!v.every(Number.isFinite)) { dropped++; continue; }
    r.push(v[0]); th.push(v[1]); ph.push(v[2]);
    if (itime >= 0) t.push(Number(c[itime]));
    if (ipass >= 0) pass.push(c[ipass]);
    if (hasM) { mbr.push(Number(c[im[0]])); mbt.push(Number(c[im[1]])); mbp.push(Number(c[im[2]])); }
  }
  if (!r.length) throw new Error('No valid numeric rows found.');
  const hasT = itime >= 0 && t.every(Number.isFinite);
  const measured = hasM && [mbr, mbt, mbp].every((a) => a.every(Number.isFinite))
    ? { br: mbr, bt: mbt, bp: mbp } : null;
  return { r, th, ph, t: hasT ? t : null, tLabel: hasT ? header[itime] : null, dropped, measured,
    pass: ipass >= 0 ? pass : null };
}

// The one built-in example: a Juno-like, near-polar flyby with closest approach
// 0.1 planetary radii above the surface at about 70 N (nearly parabolic,
// inclination 85 degrees), in degrees. Distances are in planetary radii, so it
// works for every body.
function exampleFlyby() {
  const q = 1.1, e = 1.05, rFar = 6, n = 361;
  const p = q * (1 + e);
  const nuMax = Math.acos((p / rFar - 1) / e);
  const d2r = Math.PI / 180;
  const inc = 85 * d2r, node = 30 * d2r, argp = 70 * d2r;
  const r = [], th = [], ph = [];
  for (let k = 0; k < n; k++) {
    const nu = -nuMax + (2 * nuMax * k) / (n - 1);
    const rr = p / (1 + e * Math.cos(nu));
    const u = argp + nu;
    const x = rr * (Math.cos(node) * Math.cos(u) - Math.sin(node) * Math.sin(u) * Math.cos(inc));
    const y = rr * (Math.sin(node) * Math.cos(u) + Math.cos(node) * Math.sin(u) * Math.cos(inc));
    const z = rr * Math.sin(u) * Math.sin(inc);
    r.push(+rr.toFixed(5));
    th.push(+(Math.acos(z / rr) / d2r).toFixed(4));
    ph.push(+(Math.atan2(y, x) / d2r).toFixed(4));
  }
  return { r, th, ph, t: null, tLabel: null, dropped: 0, measured: null, fixedUnit: 'deg', name: 'Synthetic flyby' };
}

function trajAngleUnit() {
  if (traj.fixedUnit) return traj.fixedUnit;
  if (ui['traj-angles'].value !== 'auto') return ui['traj-angles'].value;
  const [tlo, thi] = minmax(traj.th), [plo, phi] = minmax(traj.ph);
  return Math.max(-tlo, thi) > Math.PI + 1e-6 || Math.max(-plo, phi) > 2 * Math.PI + 1e-6 ? 'deg' : 'rad';
}

function prepareTrajectory() {
  traj.unit = trajAngleUnit();
  const f = traj.unit === 'deg' ? Math.PI / 180 : 1;
  const n = traj.r.length;
  traj.theta = Float64Array.from(traj.th, (v) => v * f);
  traj.phi = Float64Array.from(traj.ph, (v) => v * f);
  traj.xyz = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) {
    const r = traj.r[i], st = Math.sin(traj.theta[i]);
    traj.xyz.set([r * st * Math.cos(traj.phi[i]), r * st * Math.sin(traj.phi[i]), r * Math.cos(traj.theta[i])], 3 * i);
  }
  traj.rmaxPath = minmax(traj.r)[1];
  // The path itself is drawn by colorTrajectory once the field is known.
  const size = 0.03 * Math.max(1, Math.min(traj.rmaxPath, 12) / 3);
  trajStart.scale.setScalar(size);
  trajStart.position.fromArray(traj.xyz, 0);
  trajStart.visible = true;
  trajHover.scale.setScalar(size * 1.25);
}

// Split the trajectory into segments (separate passes): where the pass column
// changes, or where time runs backwards or jumps by far more than its usual step.
function trajSegments() {
  const n = traj.r.length;
  const cuts = [0];
  let typical = Infinity;
  if (traj.t && n > 2) {
    const dts = [];
    for (let i = 1; i < n; i++) dts.push(Math.abs(traj.t[i] - traj.t[i - 1]));
    dts.sort((a, b) => a - b);
    typical = dts[Math.floor(dts.length / 2)] || Infinity;
  }
  for (let i = 1; i < n; i++) {
    const passChange = traj.pass && traj.pass[i] !== traj.pass[i - 1];
    const timeBreak = traj.t && (traj.t[i] < traj.t[i - 1] || traj.t[i] - traj.t[i - 1] > 20 * typical);
    if (passChange || timeBreak) cuts.push(i);
  }
  cuts.push(n);
  const segs = [];
  for (let k = 0; k + 1 < cuts.length; k++) {
    segs.push({ start: cuts[k], end: cuts[k + 1], label: traj.passLabels ? traj.passLabels[k] : null });
  }
  return segs;
}

// (Re)draw the path, coloured by Br / max|Br| along it with the current
// colormap (blue negative, red positive for RdBu; follows Invert).
function colorTrajectory() {
  for (const o of trajGroup.children) o.geometry.dispose();
  trajGroup.clear();
  if (!traj || !trajB || traj.r.length < 2) return;
  const n = traj.r.length;
  let mx = 0;
  const measuredBr = traj.measured ? traj.measured.br : null;
  const brPath = measuredBr || trajB.br; // measured Br when the file has it
  for (const v of brPath) mx = Math.max(mx, Math.abs(v));
  const table = lut(ui.cmap.value, ui.invert.checked);
  const linear = new Float32Array(256 * 3);
  const c = new THREE.Color();
  for (let k = 0; k < 256; k++) {
    c.setRGB(table[3 * k] / 255, table[3 * k + 1] / 255, table[3 * k + 2] / 255, THREE.SRGBColorSpace);
    linear.set([c.r, c.g, c.b], 3 * k);
  }
  const colors = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) {
    const u = mx > 0 ? brPath[i] / mx : 0; // in [-1, 1]
    const k = Math.round(((u + 1) / 2) * 255);
    colors.set(linear.subarray(3 * k, 3 * k + 3), 3 * i);
  }
  for (const sg of trajSegments()) {
    if (sg.end - sg.start < 2) continue;
    const pos = traj.xyz.subarray(3 * sg.start, 3 * sg.end);
    const outlineGeo = new LineGeometry();
    outlineGeo.setPositions(pos);
    const outline = new Line2(outlineGeo, trajOutlineMat);
    outline.renderOrder = 1;
    const geo = new LineGeometry();
    geo.setPositions(pos);
    geo.setColors(colors.subarray(3 * sg.start, 3 * sg.end));
    const line = new Line2(geo, trajMat);
    line.renderOrder = 2;
    trajGroup.add(outline, line);
  }

  // Key in the side panel.
  const stops = [];
  for (let k = 0; k <= 8; k++) {
    const i = Math.round((k / 8) * 255);
    stops.push(`rgb(${table[3 * i]},${table[3 * i + 1]},${table[3 * i + 2]}) ${(k / 8) * 100}%`);
  }
  $('traj-key-grad').style.background = `linear-gradient(to right, ${stops.join(',')})`;
  const [fac, label] = val.unit;
  $('traj-key-text').innerHTML = `${measuredBr ? 'measured ' : ''}B<sub>r</sub> / max|B<sub>r</sub>| along the path ` +
    `(max |B<sub>r</sub>| = ${fmt(mx * fac)} ${label})`;
  requestRender();
}

function loadTrajectory(parsed, source) {
  traj = { ...parsed, source };
  prepareTrajectory();
  schedule(updateTrajectory);
  // Re-frame once the side panel has appeared (next frame, after schedule).
  requestAnimationFrame(() => { resize(); fitView(); });
}

function clearTrajectory() {
  traj = trajB = null;
  for (const o of trajGroup.children) o.geometry.dispose();
  trajGroup.clear();
  trajStart.visible = trajHover.visible = false;
  schedule(updateTrajectory);
  requestAnimationFrame(resize);
}

function setTrajInfo(text, warn = false) {
  const info = $('traj-info');
  info.hidden = !text;
  info.textContent = text || '';
  info.classList.toggle('warn', warn);
}

function updateTrajectory() {
  const has = !!traj;
  $('traj-clear').hidden = !has;
  $('traj-download').hidden = !has;
  if (!has) { setTrajInfo(''); return; }
  const n = traj.r.length;
  const br = new Float64Array(n), bt = new Float64Array(n), bp = new Float64Array(n), b = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const v = model.evalB(traj.r[i], traj.theta[i], traj.phi[i]);
    br[i] = v[0]; bt[i] = v[1]; bp[i] = v[2]; b[i] = Math.hypot(v[0], v[1], v[2]);
  }
  trajB = { br, bt, bp, b };
  colorTrajectory();

  const [fac, label] = val.unit;
  const [rlo, rhi] = minmax(traj.r);
  const [tlo, thi] = minmax(traj.theta);
  const notes = [];
  if (traj.dropped) notes.push(`${traj.dropped} invalid row${traj.dropped > 1 ? 's' : ''} skipped`);
  if (tlo < -1e-9 || thi > Math.PI + 1e-9) notes.push('theta outside 0–180°: is it latitude instead of co-latitude?');
  const inside = traj.r.filter((r) => r < 1).length;
  if (inside) notes.push(`${inside} point${inside > 1 ? 's' : ''} below the surface`);
  if (traj.planet && traj.planet !== val.planet) notes.push(`recorded at ${cap(traj.planet)}, shown in ${cap(val.planet)}'s field`);
  const unitNote = traj.name ? traj.name
    : `angles in ${traj.unit === 'deg' ? 'degrees' : 'radians'}${ui['traj-angles'].value === 'auto' ? ' (auto-detected)' : ''}`;
  setTrajInfo(`${n} points · r ${fmt(rlo)}–${fmt(rhi)} Rₚ · ${unitNote}` + (notes.length ? ` · ${notes.join(' · ')}` : ''), notes.length > 0);

  // Show the panel before drawing so the plot can measure its width.
  $('traj-block').hidden = false;
  $('side-panel').hidden = false;
  const credit = $('traj-credit');
  credit.hidden = !traj.credit;
  if (traj.credit) {
    credit.replaceChildren(traj.credit + (traj.url ? ', ' : '.'));
    if (traj.url) {
      const a = document.createElement('a');
      a.href = traj.url; a.target = '_blank'; a.rel = 'noopener';
      a.textContent = traj.url.replace('https://doi.org/', 'doi:');
      credit.append(a);
    }
  }
  const x = traj.t ? Float64Array.from(traj.t) : Float64Array.from({ length: n }, (_, i) => i);
  const scaled = { br, bt, bp, b };
  let measured = null;
  const misfit = $('traj-misfit');
  misfit.hidden = !traj.measured;
  if (traj.measured) {
    const M = traj.measured;
    measured = { br: M.br, bt: M.bt, bp: M.bp, b: M.br.map((v, i) => Math.hypot(v, M.bt[i], M.bp[i])) };
    // Optionally only near the planet, where the internal field dominates.
    const rlim = traj.misfitRmax || Infinity;
    let se = 0, sb = 0, nm = 0;
    for (let i = 0; i < n; i++) {
      if (traj.r[i] > rlim) continue;
      se += (br[i] - M.br[i]) ** 2 + (bt[i] - M.bt[i]) ** 2 + (bp[i] - M.bp[i]) ** 2;
      sb += measured.b[i] ** 2;
      nm++;
    }
    const where = Number.isFinite(rlim) ? ` within r ≤ ${+rlim.toPrecision(4)} Rₚ` : '';
    misfit.textContent = nm
      ? `RMS model − measured${where}: ${fmt(Math.sqrt(se / nm) * fac)} ${label} ` +
        `(RMS measured |B|: ${fmt(Math.sqrt(sb / nm) * fac)} ${label})`
      : `No points within r ≤ ${+rlim.toPrecision(4)} Rₚ for the misfit.`;
  }
  drawTrajectoryPlot($('traj-plot'), $('traj-legend'), $('traj-tip'), {
    x,
    xLabel: traj.t ? traj.tLabel : 'Point along trajectory',
    segments: trajSegments(),
    unit: label,
    series: TRAJ_SERIES.map((s) => ({
      ...s,
      values: Float64Array.from(scaled[s.key], (v) => v * fac),
      measured: measured ? Float64Array.from(measured[s.key], (v) => v * fac) : null,
    })),
    onHover: (i) => {
      trajHover.visible = i !== null;
      if (i !== null) trajHover.position.fromArray(traj.xyz, 3 * i);
      requestRender();
    },
  });
}

function downloadTrajectory() {
  if (!traj || !trajB) return;
  const [fac] = val.unit;
  const u = { nT: 'nT', muT: 'muT', Gauss: 'G' }[ui.units.value];
  const angle = traj.unit === 'deg' ? 'deg' : 'rad';
  const head = ['r', `theta_${angle}`, `phi_${angle}`];
  if (traj.t) head.push(traj.tLabel);
  head.push(`Br_${u}`, `Btheta_${u}`, `Bphi_${u}`, `B_${u}`);
  if (traj.measured) head.push(`Br_measured_${u}`, `Btheta_measured_${u}`, `Bphi_measured_${u}`);
  const rows = [
    `# planetMagFields: ${val.planet} ${val.model}${val.year !== null ? ` year ${val.year}` : ''}` +
      `${filterActive() ? ' (filtered)' : ''}; r in planetary radii, theta co-latitude, phi east longitude`,
    head.join(','),
  ];
  const g = (v) => +v.toPrecision(8);
  for (let i = 0; i < traj.r.length; i++) {
    const row = [traj.r[i], traj.th[i], traj.ph[i]];
    if (traj.t) row.push(traj.t[i]);
    row.push(g(trajB.br[i] * fac), g(trajB.bt[i] * fac), g(trajB.bp[i] * fac), g(trajB.b[i] * fac));
    if (traj.measured) {
      const M = traj.measured;
      row.push(g(M.br[i] * fac), g(M.bt[i] * fac), g(M.bp[i] * fac));
    }
    rows.push(row.join(','));
  }
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([rows.join('\n') + '\n'], { type: 'text/csv' }));
  a.download = `trajectory_${val.planet}_${val.model}.csv`;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function readTrajectoryText(text, source) {
  try {
    loadTrajectory(parseTrajectoryCSV(text), source);
  } catch (err) {
    setTrajInfo(err.message, true);
  }
}
$('traj-file').addEventListener('change', (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => readTrajectoryText(String(reader.result), 'file');
  reader.onerror = () => setTrajInfo('Could not read the file.', true);
  reader.readAsText(file);
  e.target.value = ''; // allow re-loading the same file
});
$('traj-apply').addEventListener('click', () => readTrajectoryText($('traj-text').value, 'paste'));
// Examples: real spacecraft data bundled by build_data.py, plus the synthetic
// flyby (which works for any planet).
const EXAMPLES = self.PMF_EXAMPLES || [];
for (const ex of EXAMPLES) $('traj-example').add(new Option(ex.name, ex.id));
$('traj-example').add(new Option('Synthetic flyby (any planet)', 'synthetic'));
function loadExample(id) {
  if (id === 'synthetic') { loadTrajectory(exampleFlyby(), 'synthetic'); return; }
  const ex = EXAMPLES.find((e) => e.id === id);
  if (!ex) return;
  if (ex.planet && val.planet !== ex.planet) {
    ui.planet.value = ex.planet;
    ui.planet.dispatchEvent(new Event('change'));
  }
  if (ex.units && UNITS[ex.units] && ui.units.value !== ex.units) {
    ui.units.value = ex.units;
    schedule(recolor);
  }
  const parsed = parseTrajectoryCSV(ex.csv);
  loadTrajectory({ ...parsed, fixedUnit: ex.angles, tLabel: ex.tLabel || parsed.tLabel, name: ex.name,
    credit: ex.credit, url: ex.url, planet: ex.planet, passLabels: ex.passes, misfitRmax: ex.misfitRmax }, id);
}
$('traj-example').addEventListener('change', (e) => {
  loadExample(e.target.value);
  e.target.selectedIndex = 0; // back to "Examples…" so the same one can be re-picked
});
$('traj-clear').addEventListener('click', clearTrajectory);
$('traj-download').addEventListener('click', downloadTrajectory);
$('traj-angles').addEventListener('change', () => { if (traj) { prepareTrajectory(); schedule(updateTrajectory); } });
// Redraw the plot at the new size when the side panel is resized.
let trajPlotWidth = 0;
new ResizeObserver(([entry]) => {
  const w = Math.round(entry.contentRect.width);
  if (traj && w && Math.abs(w - trajPlotWidth) > 2) { trajPlotWidth = w; schedule(updateTrajectory); }
}).observe(document.querySelector('.traj-plot-wrap'));

function updateMeta() {
  const m = DATA.planets[val.planet].models[val.model];
  const earth = val.planet === 'earth';
  $('year-row').hidden = !earth;
  $('coast-row').hidden = !earth;
  coastGroup.visible = earth && ui['show-coast'].checked;
  coastMat.color.set(ui['coast-color'].value);
  $('year-out').textContent = ui.year.value;
  $('r-out').textContent = val.r.toFixed(2);
  $('nlines-out').textContent = ui.nlines.value;
  $('rmax-out').textContent = val.rmax.toFixed(1);
  $('lw-out').textContent = ui.lw.value;
  $('lmin-out').textContent = ui.lmin.value;
  $('lmaxf-out').textContent = ui.lmaxf.value;
  $('mmin-out').textContent = ui.mmin.value;
  $('mmaxf-out').textContent = ui.mmaxf.value;
  $('filter-reset').hidden = !filterActive();
  if (!earth && playing) stopAnimation();
  $('lines-opts').hidden = !ui['show-lines'].checked;
  $('range-row').hidden = ui.sym.checked;
  $('sat-row').hidden = !ui.sym.checked;
  $('spectrum-block').hidden = !ui['show-spectrum'].checked;
  $('traj-block').hidden = !traj;
  $('side-panel').hidden = !ui['show-spectrum'].checked && !traj;
  axesGroup.visible = ui['show-axes'].checked;

  const ref = $('model-ref');
  ref.hidden = !m.url;
  if (m.url) ref.href = m.url;

  $('m-planet').textContent = cap(val.planet);
  const mm = $('m-model');
  mm.replaceChildren();
  if (m.url) {
    const a = document.createElement('a');
    a.href = m.url; a.target = '_blank'; a.rel = 'noopener'; a.textContent = m.name;
    mm.append(a);
  } else mm.textContent = m.name;
  $('m-lmax').textContent = model.lmax;
  $('m-tilt').textContent = `${fullModel.dipoleTilt().toFixed(2)}°`;
  $('m-r').textContent = val.r.toFixed(2);
  $('m-year-box').hidden = !earth;
  $('m-year').textContent = earth ? val.year.toFixed(1) : '';
  if (!playing) syncURL(); // avoid hammering history.replaceState every frame
}

const status = $('status');
function setStatus(text) {
  status.hidden = !text;
  status.textContent = text || '';
}

// Coalesce rapid slider input into one update per animation frame.
let firstFit = true;
const pending = new Set();
let scheduled = false;
function schedule(...steps) {
  steps.forEach((s) => pending.add(s));
  if (scheduled) return;
  scheduled = true;
  requestAnimationFrame(() => {
    scheduled = false;
    const run = (s) => pending.has(s) && s();
    run(buildModel);
    if (!pending.has(buildModel)) run(applyFilter);
    run(updateSurface);
    if (!pending.has(updateSurface)) run(recolor);
    run(updateLines);
    run(buildLines);
    run(updateSpectrum);
    run(updateTrajectory);
    updateMeta();
    pending.clear();
    if (firstFit) {
      // Frame the scene once the spectrum panel and colorbar have their real size.
      firstFit = false;
      resize();
      fitView();
    }
  });
}

const ALL = [buildModel, updateSurface, updateLines, updateSpectrum, updateTrajectory];
ui.planet.addEventListener('change', () => { fillModels(); resetFilter(); schedule(...ALL); fitView(); });
ui.model.addEventListener('change', () => { resetFilter(); schedule(...ALL); });
ui.year.addEventListener('input', () => { stopAnimation(); schedule(...ALL); });

// Filter sliders: keep min <= max by pushing the partner slider along.
const FILTER_STEPS = [applyFilter, updateSurface, updateLines, updateSpectrum, updateTrajectory];
for (const [lo, hi] of [['lmin', 'lmaxf'], ['mmin', 'mmaxf']]) {
  ui[lo].addEventListener('input', () => {
    if (+ui[lo].value > +ui[hi].value) ui[hi].value = ui[lo].value;
    schedule(...FILTER_STEPS);
  });
  ui[hi].addEventListener('input', () => {
    if (+ui[hi].value < +ui[lo].value) ui[lo].value = ui[hi].value;
    schedule(...FILTER_STEPS);
  });
}
$('filter-reset').addEventListener('click', () => { resetFilter(); schedule(...FILTER_STEPS); });

// Earth: animate the IGRF field through time, looping 1900 -> 2030.
let playing = false;
let animYear = 0;
let lastT = 0;
const ICON_PLAY = '<svg width="10" height="11" viewBox="0 0 10 11" aria-hidden="true"><path d="M0 0l10 5.5L0 11z" fill="currentColor"/></svg> Play';
const ICON_PAUSE = '<svg width="10" height="11" viewBox="0 0 10 11" aria-hidden="true"><path d="M0 0h3.5v11H0zM6.5 0H10v11H6.5z" fill="currentColor"/></svg> Pause';
ui.play.innerHTML = ICON_PLAY;
function stopAnimation() {
  if (!playing) return;
  playing = false;
  ui.play.innerHTML = ICON_PLAY;
  ui.play.classList.remove('playing');
  ui.play.setAttribute('aria-label', 'Play animation');
  syncURL();
  schedule(recolor); // back to the current year's own ±max
}
function animate(t) {
  if (!playing) return;
  const dt = lastT ? Math.min(0.1, (t - lastT) / 1000) : 0;
  lastT = t;
  animYear += dt * parseFloat(ui.speed.value);
  if (animYear > +ui.year.max) animYear = +ui.year.min;
  ui.year.value = animYear.toFixed(1);
  schedule(...ALL);
  requestAnimationFrame(animate);
}
ui.play.addEventListener('click', () => {
  if (playing) return stopAnimation();
  playing = true;
  animYear = parseFloat(ui.year.value);
  if (animYear >= +ui.year.max) animYear = +ui.year.min;
  lastT = 0;
  ui.play.innerHTML = ICON_PAUSE;
  ui.play.classList.add('playing');
  ui.play.setAttribute('aria-label', 'Pause animation');
  setStatus(null);
  requestAnimationFrame(animate);
});
ui.r.addEventListener('input', () => schedule(updateSurface, updateLines, updateSpectrum));
ui.resolution.addEventListener('change', () => schedule(updateSurface));
ui.units.addEventListener('change', () => schedule(recolor, updateTrajectory));
for (const id of ['cmap', 'invert', 'sym', 'sat', 'vmin', 'vmax']) ui[id].addEventListener('input', () => schedule(recolor));
ui['show-spectrum'].addEventListener('change', () => {
  schedule(updateSpectrum);
  requestAnimationFrame(() => { resize(); fitView(); });
});
for (const id of ['show-axes', 'show-coast', 'coast-color']) ui[id].addEventListener('input', () => { schedule(); requestRender(); });
for (const id of ['show-lines', 'nlines', 'seeding']) ui[id].addEventListener('input', () => schedule(updateLines));
ui.rmax.addEventListener('input', () => { schedule(updateLines); buildAxes(axisLimit()); });
for (const id of ['line-color', 'lw']) ui[id].addEventListener('input', () => schedule(buildLines));
document.querySelector('[data-reset="r"]').addEventListener('click', () => {
  ui.r.value = 1;
  schedule(updateSurface, updateLines, updateSpectrum);
});
$('reset-view').addEventListener('click', fitView);
$('open-side').addEventListener('click', () => document.body.classList.add('side-open'));
$('close-side').addEventListener('click', () => document.body.classList.remove('side-open'));

// ---------------------------------------------------------------------------
// View, hover readout and rendering
// ---------------------------------------------------------------------------
const axisLimit = () => Math.max(ui['show-lines'].checked ? val.rmax : 1, val.r) + 0.5;
const trajExtent = () => (traj ? Math.min(traj.rmaxPath, 24) : 0);

// Width taken by the colorbar on the right; the view is shifted left by half
// of it so the planet is centred in the free area instead of under the bar.
const colorbarSpace = () => $('colorbar').offsetWidth + 24;

function fitView() {
  const lim = Math.max(val.r * 1.4, Math.min(axisLimit(), 6) * 0.62, trajExtent() * 1.1);
  // Fit whichever field of view is tighter, so narrow (laptop + spectrum)
  // and very wide (ultrawide) viewports both frame the scene sensibly.
  const w = Math.max(1, viewport.clientWidth - colorbarSpace()), h = Math.max(1, viewport.clientHeight);
  const vHalf = THREE.MathUtils.degToRad(camera.fov / 2);
  const hHalf = Math.atan(Math.tan(vHalf) * (w / h));
  const dist = lim / Math.sin(Math.min(vHalf, hHalf)) * 0.75;
  camera.position.set(1, 1, 0.55).normalize().multiplyScalar(dist);
  controls.target.set(0, 0, 0);
  controls.update();
  requestRender();
}

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();
const tooltip = $('tooltip');
let hoverEvent = null;
canvas.addEventListener('pointermove', (e) => {
  if (!hoverEvent) requestAnimationFrame(() => { hover(hoverEvent); hoverEvent = null; });
  hoverEvent = e;
});
canvas.addEventListener('pointerleave', () => { tooltip.hidden = true; });

function hover(e) {
  if (!e || e.buttons) { tooltip.hidden = true; return; }
  const rect = canvas.getBoundingClientRect();
  pointer.set(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
  raycaster.setFromCamera(pointer, camera);
  const hit = raycaster.intersectObject(surface)[0];
  if (!hit) { tooltip.hidden = true; return; }
  const p = hit.point;
  const r = p.length();
  const th = Math.acos(p.z / r), ph = Math.atan2(p.y, p.x);
  const [fac, label] = val.unit;
  const br = model.evalB(r, th, ph)[0] * fac;
  const lat = 90 - THREE.MathUtils.radToDeg(th);
  const lon = THREE.MathUtils.radToDeg(ph);
  tooltip.innerHTML =
    `B<sub>r</sub>: <b>${fmt(br)} ${label}</b><br>` +
    `Lat ${Math.abs(lat).toFixed(1)}°${lat >= 0 ? 'N' : 'S'} · Lon ${Math.abs(lon).toFixed(1)}°${lon >= 0 ? 'E' : 'W'}`;
  tooltip.hidden = false;
  const x = e.clientX - rect.left + 14, y = e.clientY - rect.top + 14;
  tooltip.style.left = `${Math.min(x, rect.width - tooltip.offsetWidth - 8)}px`;
  tooltip.style.top = `${Math.min(y, rect.height - tooltip.offsetHeight - 8)}px`;
}

function resize() {
  const w = viewport.clientWidth, h = viewport.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.setViewOffset(w, h, colorbarSpace() / 2, 0, w, h);
  camera.updateProjectionMatrix();
  lineMat.resolution.set(w, h);
  trajMat.resolution.set(w, h);
  trajOutlineMat.resolution.set(w, h);
  requestRender();
}
new ResizeObserver(resize).observe(viewport);

let renderQueued = false;
function requestRender() {
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    if (controls.update()) requestRender(); // keep going while damping settles
    renderer.render(scene, camera);
  });
}

// Initial draw.
ui['show-spectrum'].checked = params.get('spectrum') === '1';
buildAxes(axisLimit());
resize();
schedule(...ALL);
if (params.has('traj')) loadExample(params.get('traj') === 'example' ? 'synthetic' : params.get('traj'));
})().catch((err) => {
  console.error(err);
  if (err instanceof LoadError) return showLoadError(err.title, err.html);
  const offline = /import|fetch|module/i.test(String(err && err.message));
  showLoadError('The app could not start',
    offline
      ? '<p>three.js could not be loaded from the CDN (cdn.jsdelivr.net). Check your internet connection and reload.</p>'
      : `<p>${String(err && err.message).replace(/</g, '&lt;')}</p>`);
});

function LoadError(title, html) { this.title = title; this.html = html; }
