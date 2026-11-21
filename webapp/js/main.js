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
const { FieldModel, seedPoints, fibonacciSeeds, traceLine, COLORMAPS, PLASMA, lut, sample, drawSpectrum } = self.PMF;

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
  'show-coast', 'coast-color', 'cmap', 'invert', 'sym', 'vmin', 'vmax', 'show-lines', 'nlines',
  'seeding', 'rmax', 'line-color', 'lw', 'lmin', 'lmaxf', 'mmin', 'mmaxf', 'play', 'speed']) ui[id] = $(id);

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
  let absmax = 0;
  for (const v of brGrid) absmax = Math.max(absmax, Math.abs(v));
  absmax *= fac;
  if (ui.sym.checked) {
    // While animating, use the largest |Br| over the whole time span so the
    // colours mean the same thing in every frame.
    if (playing) absmax = animAbsMax() * fac;
    range = [-absmax, absmax];
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
  drawColorbar(table, label);
  requestRender();
}

function fmt(v) {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e4 || a < 1e-2)) return v.toExponential(2);
  return v.toFixed(a >= 100 ? 0 : a >= 10 ? 1 : 2);
}

function drawColorbar(table, label) {
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
    d.textContent = fmt(range[0] + (k / 4) * (range[1] - range[0]));
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
  $('spectrum-panel').hidden = !ui['show-spectrum'].checked;
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

const ALL = [buildModel, updateSurface, updateLines, updateSpectrum];
ui.planet.addEventListener('change', () => { fillModels(); resetFilter(); schedule(...ALL); fitView(); });
ui.model.addEventListener('change', () => { resetFilter(); schedule(...ALL); });
ui.year.addEventListener('input', () => { stopAnimation(); schedule(...ALL); });

// Filter sliders: keep min <= max by pushing the partner slider along.
const FILTER_STEPS = [applyFilter, updateSurface, updateLines, updateSpectrum];
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
ui.units.addEventListener('change', () => schedule(recolor));
for (const id of ['cmap', 'invert', 'sym', 'vmin', 'vmax']) ui[id].addEventListener('input', () => schedule(recolor));
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

// Width taken by the colorbar on the right; the view is shifted left by half
// of it so the planet is centred in the free area instead of under the bar.
const colorbarSpace = () => $('colorbar').offsetWidth + 24;

function fitView() {
  const lim = Math.max(val.r * 1.4, Math.min(axisLimit(), 6) * 0.62);
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
