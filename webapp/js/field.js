// Classic script (no ES modules) so the app also runs from file://.
// Exposes self.PMF.{FieldModel, traceLine, ...}. The module function's source is
// also exposed (PMF.fieldSource) so the page can start a Blob worker with it,
// which works even on file:// pages where loading a worker file is blocked.
function pmfFieldModule(root) {
  'use strict';

  // Spherical-harmonic evaluation of planetary magnetic fields.
  //
  // Browser port of planetmagfields.libgauss. Coefficients use Schmidt
  // semi-normalised Legendre functions without the Condon-Shortley phase, with
  // longitude measured eastward from the prime meridian. Lengths are in planetary
  // radii and fields in nT.
  //
  // Internal potential:  V = a * sum_l (a/r)^(l+1) sum_m (g cos(m phi) + h sin(m phi)) P_l^m(cos theta)
  //   Br     = sum_l (l+1) r^-(l+2) sum_m (g cos + h sin) P
  //   Btheta = -sum_l r^-(l+2) sum_m (g cos + h sin) dP/dtheta
  //   Bphi   = sum_l r^-(l+2) sum_m m (g sin - h cos) P / sin(theta)

  const idx = (l, m) => (l * (l + 1)) / 2 + m;

  /**
   * Schmidt semi-normalised P_l^m(cos theta) and dP_l^m/dtheta for all
   * 0 <= m <= l <= lmax, written into the triangular arrays P and dP.
   */
  function legendre(lmax, theta, P, dP) {
    const c = Math.cos(theta);
    const s = Math.sin(theta);
    P[0] = 1;
    dP[0] = 0;
    for (let m = 0; m <= lmax; m++) {
      const imm = idx(m, m);
      if (m > 0) {
        const ip = idx(m - 1, m - 1);
        // P_1^1 = sin(theta); P_m^m = sqrt((2m-1)/(2m)) sin(theta) P_{m-1}^{m-1}
        const f = m === 1 ? 1 : Math.sqrt((2 * m - 1) / (2 * m));
        P[imm] = f * s * P[ip];
        dP[imm] = f * (c * P[ip] + s * dP[ip]);
      }
      for (let l = m + 1; l <= lmax; l++) {
        const i1 = idx(l - 1, m);
        const k = Math.sqrt(l * l - m * m);
        const il = idx(l, m);
        if (l === m + 1) {
          P[il] = ((2 * l - 1) * c * P[i1]) / k;
          dP[il] = ((2 * l - 1) * (c * dP[i1] - s * P[i1])) / k;
        } else {
          const i2 = idx(l - 2, m);
          const k2 = Math.sqrt((l - 1) * (l - 1) - m * m);
          P[il] = ((2 * l - 1) * c * P[i1] - k2 * P[i2]) / k;
          dP[il] = ((2 * l - 1) * (c * dP[i1] - s * P[i1]) - k2 * dP[i2]) / k;
        }
      }
    }
  }

  class FieldModel {
    /**
     * @param {number} lmax  maximum degree
     * @param {ArrayLike<number>} g  triangular g_l^m (nT), index idx(l,m)
     * @param {ArrayLike<number>} h  triangular h_l^m (nT)
     * @param {number} mmax  0 for axisymmetric models
     */
    constructor(lmax, g, h, mmax = lmax) {
      this.lmax = lmax;
      this.mmax = mmax;
      const n = idx(lmax, lmax) + 1;
      this.g = new Float64Array(n);
      this.h = new Float64Array(n);
      for (let l = 1; l <= lmax; l++) {
        for (let m = 0; m <= (mmax > 0 ? l : 0); m++) {
          const i = idx(l, m);
          this.g[i] = g[i] || 0;
          this.h[i] = h[i] || 0;
        }
      }
      this._P = new Float64Array(n);
      this._dP = new Float64Array(n);
      this._cm = new Float64Array(lmax + 1);
      this._sm = new Float64Array(lmax + 1);
    }

    /** IGRF-style time-dependent model: linear interpolation/extrapolation between epochs. */
    static fromEpochs(lmax, years, gs, hs, year, mmax = lmax) {
      let k = 0;
      while (k < years.length - 2 && year > years[k + 1]) k++;
      const t = (year - years[k]) / (years[k + 1] - years[k]);
      const lerp = (a, b) => a.map((v, i) => v + t * (b[i] - v));
      return new FieldModel(lmax, lerp(gs[k], gs[k + 1]), lerp(hs[k], hs[k + 1]), mmax);
    }

    /**
     * Copy keeping only degrees lmin..lmax and orders mmin..mmax, like
     * libgauss.filt_Gauss with lCutMin/lCutMax/mmin/mmax.
     */
    filtered(lmin, lmax, mmin, mmax) {
      const f = new FieldModel(this.lmax, this.g, this.h, this.mmax);
      for (let l = 1; l <= this.lmax; l++) {
        for (let m = 0; m <= l; m++) {
          if (l < lmin || l > lmax || (this.mmax > 0 && (m < mmin || m > mmax))) {
            f.g[idx(l, m)] = 0;
            f.h[idx(l, m)] = 0;
          }
        }
      }
      return f;
    }

    /** [Br, Btheta, Bphi] in nT at (r, theta, phi). */
    evalB(r, theta, phi) {
      const { lmax, g, h, _P: P, _dP: dP, _cm: cm, _sm: sm } = this;
      legendre(lmax, theta, P, dP);
      for (let m = 0; m <= lmax; m++) {
        cm[m] = Math.cos(m * phi);
        sm[m] = Math.sin(m * phi);
      }
      const s = Math.max(Math.sin(theta), 1e-10);
      let br = 0, bt = 0, bp = 0;
      let rl = 1 / (r * r); // r^-(l+2) for l = 0
      for (let l = 1; l <= lmax; l++) {
        rl /= r;
        let a = 0, b = 0, c = 0;
        const mtop = this.mmax > 0 ? l : 0;
        for (let m = 0; m <= mtop; m++) {
          const i = idx(l, m);
          const gc = g[i] * cm[m] + h[i] * sm[m];
          a += gc * P[i];
          b += gc * dP[i];
          if (m > 0) c += m * (g[i] * sm[m] - h[i] * cm[m]) * P[i];
        }
        br += (l + 1) * rl * a;
        bt -= rl * b;
        bp += (rl * c) / s;
      }
      return [br, bt, bp];
    }

    /** Cartesian field [Bx, By, Bz] at Cartesian point (x, y, z). */
    evalBxyz(x, y, z) {
      const rho = Math.hypot(x, y);
      const r = Math.hypot(rho, z);
      const theta = Math.atan2(rho, z);
      const phi = Math.atan2(y, x);
      const [br, bt, bp] = this.evalB(r, theta, phi);
      const st = Math.sin(theta), ct = Math.cos(theta);
      const sp = Math.sin(phi), cp = Math.cos(phi);
      return [
        br * st * cp + bt * ct * cp - bp * sp,
        br * st * sp + bt * ct * sp + bp * cp,
        br * ct - bt * st,
      ];
    }

    /**
     * Br on a cell-centred (theta, phi) grid: row j at theta = (j+0.5) pi / ntheta,
     * column i at phi = (i+0.5) 2pi / nphi. Returns Float32Array of ntheta*nphi.
     */
    gridBr(r, nphi, ntheta) {
      const { lmax, g, h, _P: P, _dP: dP } = this;
      const mtopMax = this.mmax > 0 ? lmax : 0;
      const cosm = new Float64Array((mtopMax + 1) * nphi);
      const sinm = new Float64Array((mtopMax + 1) * nphi);
      for (let i = 0; i < nphi; i++) {
        const phi = ((i + 0.5) * 2 * Math.PI) / nphi;
        for (let m = 0; m <= mtopMax; m++) {
          cosm[m * nphi + i] = Math.cos(m * phi);
          sinm[m * nphi + i] = Math.sin(m * phi);
        }
      }
      const rl = new Float64Array(lmax + 1);
      for (let l = 1; l <= lmax; l++) rl[l] = (l + 1) * Math.pow(r, -(l + 2));
      const A = new Float64Array(mtopMax + 1);
      const B = new Float64Array(mtopMax + 1);
      const out = new Float32Array(nphi * ntheta);
      for (let j = 0; j < ntheta; j++) {
        legendre(lmax, ((j + 0.5) * Math.PI) / ntheta, P, dP);
        A.fill(0);
        B.fill(0);
        for (let l = 1; l <= lmax; l++) {
          for (let m = 0; m <= Math.min(l, mtopMax); m++) {
            const i = idx(l, m);
            A[m] += rl[l] * g[i] * P[i];
            B[m] += rl[l] * h[i] * P[i];
          }
        }
        for (let i = 0; i < nphi; i++) {
          let v = 0;
          for (let m = 0; m <= mtopMax; m++) v += A[m] * cosm[m * nphi + i] + B[m] * sinm[m * nphi + i];
          out[j * nphi + i] = v;
        }
      }
      return out;
    }

    /** Lowes spectrum and energy partition at radius r (port of libgauss.get_spec). */
    spectrum(r = 1) {
      const { lmax, g, h } = this;
      const E = new Float64Array(lmax + 1);
      let symm = 0, antisymm = 0, axi = 0;
      for (let l = 1; l <= lmax; l++) {
        const f = (l + 1) * Math.pow(r, -2 * l - 4);
        const mtop = this.mmax > 0 ? l : 0;
        for (let m = 0; m <= mtop; m++) {
          const i = idx(l, m);
          const e = f * (g[i] * g[i] + h[i] * h[i]);
          E[l] += e;
          if (this.mmax === 0) continue;
          if ((l - m) % 2 === 0) symm += e; else antisymm += e;
        }
        axi += f * (g[idx(l, 0)] ** 2 + h[idx(l, 0)] ** 2);
      }
      let tot = 0;
      for (let l = 1; l <= lmax; l++) tot += E[l];
      if (this.mmax === 0) {
        // Matches get_spec for axisymmetric models: sums over even / odd l.
        symm = 0; antisymm = 0;
        for (let l = 0; l <= lmax; l++) (l % 2 === 0 ? (symm += E[l]) : (antisymm += E[l]));
      }
      const e10 = this.mmax === 0 ? E[1] : 2 * Math.pow(r, -6) * g[idx(1, 0)] ** 2;
      return {
        E, tot, symm, antisymm, axi,
        dipolarity: e10 / tot,
        dipTot: E[1] / tot,
      };
    }

    /** Dipole tilt in degrees, as in libgauss.get_dipole_tilt. */
    dipoleTilt() {
      if (this.mmax === 0) return 0;
      const g10 = this.g[idx(1, 0)], g11 = this.g[idx(1, 1)], h11 = this.h[idx(1, 1)];
      return (Math.atan(Math.hypot(g11, h11) / g10) * 180) / Math.PI;
    }
  }

  /** Small seeded PRNG so seeding is reproducible across sessions. */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return () => {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /**
   * Fibonacci (golden-spiral) lattice: n near-uniformly spaced [theta, phi]
   * seeds on the sphere, deterministic and free of clustering.
   */
  function fibonacciSeeds(n) {
    const golden = Math.PI * (3 - Math.sqrt(5));
    const seeds = [];
    for (let k = 0; k < n; k++) {
      const z = 1 - (2 * k + 1) / n;
      seeds.push([Math.acos(z), (k * golden) % (2 * Math.PI)]);
    }
    return seeds;
  }

  /**
   * Seed points for field lines, following the Streamlit app's strategy:
   * 40% weighted by |Br| in positive regions, 40% in negative regions, the rest
   * uniform on the sphere. Weighted sampling without replacement uses
   * Efraimidis-Spirakis keys u^(1/w).
   */
  function seedPoints(br, nphi, ntheta, nLines, seed = 42) {
    const rand = mulberry32(seed);
    const nPos = Math.max(1, Math.floor(0.4 * nLines));
    const nNeg = Math.max(1, Math.floor(0.4 * nLines));
    const nUni = Math.max(0, nLines - nPos - nNeg);
    const seeds = [];
    const pick = (sign, k) => {
      const keys = [];
      for (let n = 0; n < br.length; n++) {
        const w = sign * br[n];
        if (w > 0) keys.push([Math.pow(rand(), 1 / w), n]);
      }
      keys.sort((a, b) => b[0] - a[0]);
      for (const [, n] of keys.slice(0, k)) {
        const j = Math.floor(n / nphi), i = n % nphi;
        seeds.push([((j + 0.5) * Math.PI) / ntheta, ((i + 0.5) * 2 * Math.PI) / nphi]);
      }
    };
    pick(1, nPos);
    pick(-1, nNeg);
    for (let k = 0; k < nUni; k++) {
      seeds.push([Math.acos(2 * rand() - 1), 2 * Math.PI * rand()]);
    }
    return seeds;
  }

  /**
   * Trace one field line through a seed on the sphere of radius rSurf, in both
   * directions, with fixed-step RK4 on the unit field direction (step scales
   * with distance). Stops at rSurf, rMax or maxLength. Returns a flat
   * Float32Array of xyz points, or null if the line is too short.
   */
  function traceLine(model, theta0, phi0, rSurf, rMax, maxLength = 80) {
    const r0 = rSurf + 0.02;
    const start = [
      r0 * Math.sin(theta0) * Math.cos(phi0),
      r0 * Math.sin(theta0) * Math.sin(phi0),
      r0 * Math.cos(theta0),
    ];
    const dir = (p, sgn) => {
      const b = model.evalBxyz(p[0], p[1], p[2]);
      const n = Math.hypot(b[0], b[1], b[2]) || 1;
      return [(sgn * b[0]) / n, (sgn * b[1]) / n, (sgn * b[2]) / n];
    };
    const half = (sgn) => {
      const pts = [start];
      let p = start, len = 0;
      while (len < maxLength) {
        const r = Math.hypot(p[0], p[1], p[2]);
        const ds = Math.min(0.1, Math.max(0.004, 0.015 * r));
        const k1 = dir(p, sgn);
        const k2 = dir([p[0] + 0.5 * ds * k1[0], p[1] + 0.5 * ds * k1[1], p[2] + 0.5 * ds * k1[2]], sgn);
        const k3 = dir([p[0] + 0.5 * ds * k2[0], p[1] + 0.5 * ds * k2[1], p[2] + 0.5 * ds * k2[2]], sgn);
        const k4 = dir([p[0] + ds * k3[0], p[1] + ds * k3[1], p[2] + ds * k3[2]], sgn);
        const q = [0, 1, 2].map((c) => p[c] + (ds / 6) * (k1[c] + 2 * k2[c] + 2 * k3[c] + k4[c]));
        const rq = Math.hypot(q[0], q[1], q[2]);
        len += ds;
        if (rq < rSurf || rq > rMax) {
          // Land exactly on the bounding sphere so lines meet the surface.
          const rb = rq < rSurf ? rSurf : rMax;
          const t = (rb - r) / (rq - r);
          pts.push([0, 1, 2].map((c) => p[c] + t * (q[c] - p[c])));
          break;
        }
        pts.push(q);
        p = q;
      }
      return pts;
    };
    const fwd = half(1);
    const bwd = half(-1);
    const all = bwd.reverse().concat(fwd.slice(1));
    if (all.length < 5) return null;
    const out = new Float32Array(all.length * 3);
    all.forEach((p, n) => out.set(p, 3 * n));
    return out;
  }

  root.PMF = Object.assign(root.PMF || {}, { idx, legendre, FieldModel, mulberry32, fibonacciSeeds, seedPoints, traceLine });
}
pmfFieldModule(self);
self.PMF.fieldSource = pmfFieldModule.toString();
