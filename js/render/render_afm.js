// js/render/render_afm.js
// (split from renderer.js, Wave R1)

import {
  get_covalent_radius,
  guess_bonds_by_distance,
} from "../bonds/bonds.js";
import {
  compute_scaled_coordinates,
  newFloatImage,
  clamp255,
  gaussianBlurFloat,
  getGaussian1D_cached,
  draw_gaussian_line,
  draw_bond_glow_line,
  draw_scale_bar,
} from "./render_common.js";

// AFM-like surface render (Wave 2)
//  - Height/topography map (max of Gaussian bumps)
//  - Edge enhancement (simple Laplacian)
//  - Optional single-line "bond contrast" overlay (bt ignored; legacy profile with glow; no multi-lines)
//  - NOTE: DoF is intentionally ignored in AFM mode (surface mode).
// ---------------------------

function draw_bonds_single(img, coords, bonds, atoms, opts = {}) {
  // AFM bonds: "legacy" TEM bond profile (core + glow), but ALWAYS single-line (bt ignored).
  // Polarity is handled by the caller (AFM adds this layer to brighten; global invert stays global).
  const {
    wave_width_px = 8.0,
    wave_amp = 1.0,
    focal_z = 0.0,
    hide_front = false,
    z_view = null,
  } = opts;

  const halfw = Math.max(1.0, (Number(wave_width_px) || 8.0) * 0.5);
  const amp = Math.max(0.0, Number(wave_amp) || 0.0);
  if (!(amp > 0)) return;

  for (let k = 0; k < bonds.length; k++) {
    const b = bonds[k];
    const i = b[0],
      j = b[1];
    const p1 = coords[i],
      p2 = coords[j];
    if (!p1 || !p2) continue;

    // hide_front: respect focal_z for educational depth slicing
    if (hide_front) {
      const zi = z_view ? (z_view[i] ?? 0) : (atoms[i].z ?? 0);
      const zj = z_view ? (z_view[j] ?? 0) : (atoms[j].z ?? 0);
      if (zi > focal_z + 1e-9 || zj > focal_z + 1e-9) continue;
    }

    const Z1 = atoms[i]?.Z ?? 0;
    const Z2 = atoms[j]?.Z ?? 0;
    const Zavg = 0.5 * (Z1 + Z2);

    // Legacy intensity profile (same as draw_bonds), scaled by wave_amp
    const base_int = Math.min(255.0, Math.pow(Zavg, 0.9));
    const glow_int = Math.min(255.0, Math.pow(Zavg, 1.2) * 2.0);
    const bond_intensity = Math.min(255.0, base_int * amp);
    const glow_intensity = Math.min(255.0, glow_int * amp);

    // bbox culling (single line only => no offset padding)
    const pad = halfw + 2;
    const minx = Math.min(p1[0], p2[0]) - pad;
    const maxx = Math.max(p1[0], p2[0]) + pad;
    const miny = Math.min(p1[1], p2[1]) - pad;
    const maxy = Math.max(p1[1], p2[1]) + pad;
    if (maxx < 0 || minx >= img.w || maxy < 0 || miny >= img.h) continue;

    // Always one bridge (ignore bt)
    draw_gaussian_line(img, p1, p2, bond_intensity, halfw);
    draw_bond_glow_line(
      img,
      p1,
      p2,
      glow_intensity * 0.85,
      Math.max(0.5, 0.6 * halfw),
    );
  }
}


// Element overrides helper (Z->symbol fallback)
const _EO_Z2SYM = [
  null,
  "H","He",
  "Li","Be","B","C","N","O","F","Ne",
  "Na","Mg","Al","Si","P","S","Cl","Ar",
  "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
  "In","Sn","Sb","Te","I","Xe",
  "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
  "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn",
  "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
  "Nh","Fl","Mc","Lv","Ts","Og"
];

function _EO_atomSymbol(a) {
  if (!a) return null;

  // Prefer atomic number -> symbol (SMILES/RDKit can mis-populate .sym)
  let Z = a.Z;
  if (!Number.isFinite(Z)) Z = parseInt(Z, 10);
  if (Number.isFinite(Z)) {
    const zi = Z | 0;
    if (zi > 0 && zi < _EO_Z2SYM.length) {
      const zs = _EO_Z2SYM[zi];
      if (zs) return zs;
    }
  }

  let s = a.sym || a.el || a.symbol || a.element;
  if (typeof s === "string") {
    s = s.trim();
    if (s) return s;
  }
  return null;
}

function _EO_getMul(ov, sym) {
  let sizeMul = 1.0;
  let darkMul = 1.0;
  if (ov && sym && ov[sym]) {
    const it = ov[sym];
    if (it && Number.isFinite(it.size)) sizeMul = it.size;
    if (it && Number.isFinite(it.dark)) darkMul = it.dark;
  }
  return [sizeMul, darkMul];
}

function _clamp01(x) {
  return x <= 0 ? 0 : x >= 1 ? 1 : x;
}

function _smoothstep(e0, e1, x) {
  if (!(e1 > e0)) return x >= e1 ? 1 : 0;
  const t = _clamp01((x - e0) / (e1 - e0));
  return t * t * (3 - 2 * t);
}

function _normalizePositive(src, pow = 1.0) {
  let maxv = 0.0;
  for (let i = 0; i < src.length; i++) {
    const v = src[i];
    if (v > maxv) maxv = v;
  }
  const out = new Float32Array(src.length);
  if (!(maxv > 1e-12)) return out;
  const inv = 1.0 / maxv;
  if (Math.abs(pow - 1.0) < 1e-9) {
    for (let i = 0; i < src.length; i++) {
      const v = src[i] * inv;
      out[i] = v > 0 ? v : 0;
    }
  } else {
    for (let i = 0; i < src.length; i++) {
      const v = src[i] * inv;
      out[i] = v > 0 ? Math.pow(v, pow) : 0;
    }
  }
  return out;
}


function _defaultAFMRadialProfile() {
  return {
    start_A: 3.2,
    mid_A: 4.8,
    full_A: 7.2,
    gamma_inner: 0.48,
    gamma_outer: 0.72,
    tail_mix: 0.42,
    bond_gain: 2.8,
    atom_dark_gain: 0.38,
    blur_mix: 0.24,
    blur_sigma_px: 2.5,
    base_bond_gain: 1.08,
  };
}

function _sanitizeAFMRadialProfile(p) {
  const d = _defaultAFMRadialProfile();
  const out = { ...d, ...(p || {}) };

  if (!Number.isFinite(out.start_A)) out.start_A = d.start_A;
  if (!Number.isFinite(out.mid_A)) out.mid_A = d.mid_A;
  if (!Number.isFinite(out.full_A)) out.full_A = d.full_A;
  if (!Number.isFinite(out.gamma_inner)) out.gamma_inner = d.gamma_inner;
  if (!Number.isFinite(out.gamma_outer)) out.gamma_outer = d.gamma_outer;
  if (!Number.isFinite(out.tail_mix)) out.tail_mix = d.tail_mix;
  if (!Number.isFinite(out.bond_gain)) out.bond_gain = d.bond_gain;
  if (!Number.isFinite(out.atom_dark_gain)) out.atom_dark_gain = d.atom_dark_gain;
  if (!Number.isFinite(out.blur_mix)) out.blur_mix = d.blur_mix;
  if (!Number.isFinite(out.blur_sigma_px)) out.blur_sigma_px = d.blur_sigma_px;
  if (!Number.isFinite(out.base_bond_gain)) out.base_bond_gain = d.base_bond_gain;

  out.start_A = Math.max(0.0, out.start_A);
  out.mid_A = Math.max(out.start_A + 1e-6, out.mid_A);
  out.full_A = Math.max(out.mid_A + 1e-6, out.full_A);

  out.gamma_inner = Math.max(0.15, out.gamma_inner);
  out.gamma_outer = Math.max(0.15, out.gamma_outer);

  out.tail_mix = _clamp01(out.tail_mix);
  out.bond_gain = Math.max(0.0, out.bond_gain);
  out.atom_dark_gain = Math.max(0.0, out.atom_dark_gain);
  out.blur_mix = _clamp01(out.blur_mix);
  out.blur_sigma_px = Math.max(0.0, out.blur_sigma_px);
  out.base_bond_gain = Math.max(0.0, out.base_bond_gain);

  return out;
}

function _ensureAFMConsoleControls() {
  const g = globalThis;
  if (!Number.isFinite(g.__AFM_RADIAL_GAIN)) g.__AFM_RADIAL_GAIN = 1.0;
  if (!Number.isFinite(g.__AFM_RADIAL_BOND_GAIN)) g.__AFM_RADIAL_BOND_GAIN = 0.0;

  g.__AFM_RADIAL_PROFILE = _sanitizeAFMRadialProfile(g.__AFM_RADIAL_PROFILE);

  if (typeof g.__setAFMRadial !== "function") {
    g.__setAFMRadial = function (gain = 1.0, bondGain = g.__AFM_RADIAL_BOND_GAIN) {
      const gg = Number(gain);
      const bb = Number(bondGain);
      if (Number.isFinite(gg)) g.__AFM_RADIAL_GAIN = gg;
      if (Number.isFinite(bb)) g.__AFM_RADIAL_BOND_GAIN = bb;
      return {
        radial: g.__AFM_RADIAL_GAIN,
        bond: g.__AFM_RADIAL_BOND_GAIN,
        profile: g.__AFM_RADIAL_PROFILE,
      };
    };
  }

  if (typeof g.__setAFMRadialProfile !== "function") {
    g.__setAFMRadialProfile = function (patch = {}) {
      g.__AFM_RADIAL_PROFILE = _sanitizeAFMRadialProfile({
        ...g.__AFM_RADIAL_PROFILE,
        ...(patch || {}),
      });
      return g.__AFM_RADIAL_PROFILE;
    };
  }
}

function _resolveAFMRadialControls(opts = {}) {
  _ensureAFMConsoleControls();
  const g = globalThis;
  let gain = Number.isFinite(opts?.afm_radial_gain)
    ? Number(opts.afm_radial_gain)
    : Number(g.__AFM_RADIAL_GAIN);
  let bondGain = Number.isFinite(opts?.afm_radial_bond_gain)
    ? Number(opts.afm_radial_bond_gain)
    : Number(g.__AFM_RADIAL_BOND_GAIN);
  if (!Number.isFinite(gain)) gain = 1.0;
  if (!Number.isFinite(bondGain)) bondGain = 0.0;
  if (gain < 0) gain = 0;
  if (bondGain < 0) bondGain = 0;

  const optProfile =
    opts && typeof opts.afm_radial_profile === "object" && opts.afm_radial_profile
      ? opts.afm_radial_profile
      : null;
  const profile = _sanitizeAFMRadialProfile(
    optProfile ? { ...g.__AFM_RADIAL_PROFILE, ...optProfile } : g.__AFM_RADIAL_PROFILE,
  );
  return { gain, bondGain, profile };
}

function _radialProfileValue(distA, profile) {
  if (!(distA > profile.start_A)) return 0.0;

  const s = profile.start_A;
  const m = profile.mid_A;
  const f = profile.full_A;

  const innerDen = Math.max(1e-6, m - s);
  const outerDen = Math.max(1e-6, f - m);

  // Continuous two-stage radial ramp:
  //  - stage 1 (start -> mid): early growth, but with softened curvature
  //  - stage 2 (mid -> full): continues from a high shoulder instead of jumping
  // This keeps the effect clearly visible around ~5 Å while removing the mid-point step.
  const shoulder = 1.0 - 0.25 * profile.tail_mix;

  let v;
  if (distA <= m) {
    const t = _clamp01((distA - s) / innerDen);
    const eased = _smoothstep(0.0, 1.0, t);
    v = shoulder * Math.pow(eased, profile.gamma_inner);
  } else {
    const t = _clamp01((distA - m) / outerDen);
    const eased = _smoothstep(0.0, 1.0, t);
    v =
      shoulder +
      (1.0 - shoulder) * Math.pow(eased, profile.gamma_outer);
  }
  return _clamp01(v);
}


function _afmBuildBondAdjacency(nAtoms, bonds) {
  const adj = new Array(nAtoms);
  for (let i = 0; i < nAtoms; i++) adj[i] = [];
  if (!bonds || !bonds.length) return adj;
  for (let k = 0; k < bonds.length; k++) {
    const b = bonds[k];
    if (!b || b.length < 2) continue;
    const i = b[0] | 0;
    const j = b[1] | 0;
    if (i < 0 || j < 0 || i >= nAtoms || j >= nAtoms || i === j) continue;
    adj[i].push(j);
    adj[j].push(i);
  }
  return adj;
}

function _afmAtomLinearityScores(coords, adj) {
  const out = new Float32Array(adj.length);
  for (let i = 0; i < adj.length; i++) {
    const ns = adj[i];
    if (!ns || ns.length < 2) continue;
    const pi = coords[i];
    if (!pi) continue;

    let best = 0.0;
    for (let a = 0; a < ns.length; a++) {
      const pa = coords[ns[a]];
      if (!pa) continue;
      const ax = pa[0] - pi[0];
      const ay = pa[1] - pi[1];
      const al = Math.hypot(ax, ay);
      if (!(al > 1e-6)) continue;

      for (let b = a + 1; b < ns.length; b++) {
        const pb = coords[ns[b]];
        if (!pb) continue;
        const bx = pb[0] - pi[0];
        const by = pb[1] - pi[1];
        const bl = Math.hypot(bx, by);
        if (!(bl > 1e-6)) continue;

        const dot = (ax * bx + ay * by) / (al * bl);
        const col = Math.abs(dot);
        if (col > best) best = col;
      }
    }

    if (ns.length >= 3) best = Math.max(best, 0.55);
    out[i] = best;
  }
  return out;
}

function _drawAFMGaussianSpot(arr, W, H, x0, y0, amp, sigma) {
  if (!(amp > 1e-8) || !(sigma > 1e-6)) return;

  const R = Math.ceil(4 * sigma);
  if (x0 + R < 0 || x0 - R >= W || y0 + R < 0 || y0 - R >= H) return;

  const x0i = Math.floor(x0);
  const y0i = Math.floor(y0);
  const fracX = x0 - x0i;
  const fracY = y0 - y0i;
  const wx = getGaussian1D_cached(sigma, fracX, R);
  const wy = getGaussian1D_cached(sigma, fracY, R);

  const yStart = Math.max(0, y0i - R);
  const yEnd = Math.min(H - 1, y0i + R);
  const xStart = Math.max(0, x0i - R);
  const xEnd = Math.min(W - 1, x0i + R);

  for (let y = yStart; y <= yEnd; y++) {
    const wyv = wy[y - y0i + R];
    const row = y * W;
    for (let x = xStart; x <= xEnd; x++) {
      const g = wx[x - x0i + R] * wyv;
      if (g > 1e-10) arr[row + x] += amp * g;
    }
  }
}


function _sampleFieldBilinear(arr, W, H, x, y) {
  if (!arr || !W || !H) return 0.0;
  if (x < 0 || y < 0 || x > W - 1 || y > H - 1) return 0.0;

  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const x1 = Math.min(W - 1, x0 + 1);
  const y1 = Math.min(H - 1, y0 + 1);
  const tx = x - x0;
  const ty = y - y0;

  const p00 = y0 * W + x0;
  const p10 = y0 * W + x1;
  const p01 = y1 * W + x0;
  const p11 = y1 * W + x1;

  const v00 = arr[p00] || 0.0;
  const v10 = arr[p10] || 0.0;
  const v01 = arr[p01] || 0.0;
  const v11 = arr[p11] || 0.0;

  const vx0 = v00 + (v10 - v00) * tx;
  const vx1 = v01 + (v11 - v01) * tx;
  return vx0 + (vx1 - vx0) * ty;
}

function _sampleAFMStructureSupport(ctx, W, H, x, y) {
  const shape = ctx.shapeSupport ? _sampleFieldBilinear(ctx.shapeSupport, W, H, x, y) : 0.0;
  const h = ctx.height ? _sampleFieldBilinear(ctx.height, W, H, x, y) : 0.0;
  const mol = ctx.molMask ? _sampleFieldBilinear(ctx.molMask, W, H, x, y) : 0.0;
  const merge = ctx.mergePresence ? _sampleFieldBilinear(ctx.mergePresence, W, H, x, y) : 0.0;
  return _clamp01(Math.max(shape, 0.52 * mol, 0.50 * h, 0.54 * merge));
}

function _afmViewportConfidence(W, H, x, y, margin) {
  const m = Math.max(1e-6, Number(margin) || 0.0);
  const d = Math.min(x, y, (W - 1) - x, (H - 1) - y);
  return _smoothstep(0.0, m, d);
}

function _sampleAFMStructureSupportSafe(ctx, W, H, x, y, margin, fallback) {
  const conf = _afmViewportConfidence(W, H, x, y, margin);
  const v = _sampleAFMStructureSupport(ctx, W, H, x, y);
  const fb = Number.isFinite(fallback) ? fallback : 0.0;
  return {
    value: v * conf + fb * (1.0 - conf),
    conf,
  };
}

function _afmAtomTopologyScores(coords, adj) {
  const linearity = _afmAtomLinearityScores(coords, adj);
  const tip = new Float32Array(adj.length);
  const bend = new Float32Array(adj.length);
  const branch = new Float32Array(adj.length);
  const exposure = new Float32Array(adj.length);

  for (let i = 0; i < adj.length; i++) {
    const deg = adj[i] ? adj[i].length : 0;
    const lin = Math.max(0.0, Math.min(1.0, linearity[i] || 0.0));

    if (deg <= 0) continue;
    if (deg === 1) {
      tip[i] = 1.0;
      exposure[i] = 1.0;
      continue;
    }

    if (deg === 2) {
      const kink = Math.pow(1.0 - lin, 0.85);
      bend[i] = kink;
      exposure[i] = Math.max(0.08, 0.92 * kink);
      continue;
    }

    const multi = 0.45 + 0.55 * Math.pow(1.0 - lin, 0.75);
    branch[i] = multi;
    bend[i] = Math.max(0.42, multi * 0.92);
    exposure[i] = Math.max(0.38, multi * 0.86);
  }

  return { linearity, tip, bend, branch, exposure };
}


function _afmBfsDistances(adj, start, componentSet = null) {
  const n = adj.length;
  const dist = new Int32Array(n);
  dist.fill(-1);
  const q = new Int32Array(Math.max(1, n));
  let qs = 0;
  let qe = 0;
  dist[start] = 0;
  q[qe++] = start;
  let far = start;

  while (qs < qe) {
    const v = q[qs++];
    far = v;
    const ns = adj[v] || [];
    for (let k = 0; k < ns.length; k++) {
      const u = ns[k] | 0;
      if (u < 0 || u >= n || dist[u] >= 0) continue;
      if (componentSet && !componentSet.has(u)) continue;
      dist[u] = dist[v] + 1;
      q[qe++] = u;
    }
  }

  return { dist, far };
}

function _afmAtomGraphShapeScores(adj) {
  const n = adj.length;
  const endness = new Float32Array(n);
  const branchness = new Float32Array(n);
  const components = new Int32Array(n);
  components.fill(-1);

  let cid = 0;
  for (let i = 0; i < n; i++) {
    if (components[i] >= 0) continue;
    const stack = [i];
    const nodes = [];
    components[i] = cid;

    while (stack.length) {
      const v = stack.pop();
      nodes.push(v);
      const ns = adj[v] || [];
      for (let k = 0; k < ns.length; k++) {
        const u = ns[k] | 0;
        if (u < 0 || u >= n || components[u] >= 0) continue;
        components[u] = cid;
        stack.push(u);
      }
    }

    const comp = new Set(nodes);
    const first = nodes[0];
    const bfs0 = _afmBfsDistances(adj, first, comp);
    const bfsA = _afmBfsDistances(adj, bfs0.far, comp);
    const bfsB = _afmBfsDistances(adj, bfsA.far, comp);
    const distA = bfsA.dist;
    const distB = bfsB.dist;
    const diam = Math.max(1, distB[bfsA.far] | 0);

    for (let t = 0; t < nodes.length; t++) {
      const v = nodes[t];
      const da = distA[v] < 0 ? diam : distA[v];
      const db = distB[v] < 0 ? diam : distB[v];
      const inner = Math.min(da, db);
      const endK = _clamp01(1.0 - (2.0 * inner) / Math.max(1, diam));
      const deg = adj[v] ? adj[v].length : 0;
      endness[v] = endK;
      branchness[v] = deg <= 2 ? 0.0 : _clamp01((deg - 2) / 2.0);
    }

    cid++;
  }

  return { endness, branchness, components };
}

function _afmComputePrincipalAxes(coords, atoms) {
  let sw = 0.0;
  let sx = 0.0;
  let sy = 0.0;
  for (let i = 0; i < coords.length; i++) {
    const p = coords[i];
    const a = atoms[i];
    if (!p || !a) continue;
    const Z = (a.Z | 0) || 6;
    const w = Math.max(0.25, get_covalent_radius(Z));
    sw += w;
    sx += p[0] * w;
    sy += p[1] * w;
  }

  const cx = sw > 1e-8 ? sx / sw : 0.0;
  const cy = sw > 1e-8 ? sy / sw : 0.0;

  let sxx = 0.0;
  let syy = 0.0;
  let sxy = 0.0;
  for (let i = 0; i < coords.length; i++) {
    const p = coords[i];
    const a = atoms[i];
    if (!p || !a) continue;
    const Z = (a.Z | 0) || 6;
    const w = Math.max(0.25, get_covalent_radius(Z));
    const dx = p[0] - cx;
    const dy = p[1] - cy;
    sxx += w * dx * dx;
    syy += w * dy * dy;
    sxy += w * dx * dy;
  }

  const tr = sxx + syy;
  const det = sxx * syy - sxy * sxy;
  const disc = Math.sqrt(Math.max(0.0, 0.25 * tr * tr - det));
  const l1 = 0.5 * tr + disc;

  let ax1x = 1.0;
  let ax1y = 0.0;
  if (Math.abs(sxy) > 1e-8 || Math.abs(sxx - l1) > 1e-8) {
    ax1x = sxy;
    ax1y = l1 - sxx;
    const al = Math.hypot(ax1x, ax1y);
    if (al > 1e-8) {
      ax1x /= al;
      ax1y /= al;
    }
  }
  const ax2x = -ax1y;
  const ax2y = ax1x;

  let maxProj1 = 1e-6;
  let maxProj2 = 1e-6;
  let maxDist = 1e-6;
  for (let i = 0; i < coords.length; i++) {
    const p = coords[i];
    if (!p) continue;
    const dx = p[0] - cx;
    const dy = p[1] - cy;
    const pr1 = Math.abs(dx * ax1x + dy * ax1y);
    const pr2 = Math.abs(dx * ax2x + dy * ax2y);
    const dl = Math.hypot(dx, dy);
    if (pr1 > maxProj1) maxProj1 = pr1;
    if (pr2 > maxProj2) maxProj2 = pr2;
    if (dl > maxDist) maxDist = dl;
  }

  return { cx, cy, ax1x, ax1y, ax2x, ax2y, maxProj1, maxProj2, maxDist };
}

function _sampleFieldClamped(arr, W, H, x, y) {
  if (!arr || !W || !H) return 0.0;
  const xx = x < 0 ? 0 : x > W - 1 ? W - 1 : x;
  const yy = y < 0 ? 0 : y > H - 1 ? H - 1 : y;
  return _sampleFieldBilinear(arr, W, H, xx, yy);
}

function _buildAFMCoarseShapeField(coords, atoms, adj, topo, shapeField, W, H, angstroms_per_pixel) {
  const n = coords.length;
  const graph = _afmAtomGraphShapeScores(adj);
  const axes = _afmComputePrincipalAxes(coords, atoms);
  const terminal = new Float32Array(n);
  const corner = new Float32Array(n);
  const branch = new Float32Array(n);
  const peripheral = new Float32Array(n);
  const salient = new Float32Array(n);

  const outwardPx1 = Math.max(1.2, 1.30 / Math.max(angstroms_per_pixel, 1e-6));
  const outwardPx2 = Math.max(outwardPx1 + 0.8, 2.45 / Math.max(angstroms_per_pixel, 1e-6));

  for (let i = 0; i < n; i++) {
    const p = coords[i];
    if (!p) continue;

    const dx = p[0] - axes.cx;
    const dy = p[1] - axes.cy;
    const pr1 = Math.abs(dx * axes.ax1x + dy * axes.ax1y) / Math.max(1e-6, axes.maxProj1);
    const pr2 = Math.abs(dx * axes.ax2x + dy * axes.ax2y) / Math.max(1e-6, axes.maxProj2);
    const radial = Math.hypot(dx, dy) / Math.max(1e-6, axes.maxDist);
    const longDominance = _clamp01(1.0 - pr2 / Math.max(0.15, pr1 + 0.08));
    const graphEnd = graph.endness[i] || 0.0;
    const graphBranch = graph.branchness[i] || 0.0;
    const bend = topo.bend[i] || 0.0;
    const topoBranch = topo.branch[i] || 0.0;

    let ux = dx;
    let uy = dy;
    let ul = Math.hypot(ux, uy);
    if (!(ul > 1e-6)) {
      const sign = (dx * axes.ax1x + dy * axes.ax1y) >= 0 ? 1.0 : -1.0;
      ux = sign * axes.ax1x;
      uy = sign * axes.ax1y;
      ul = 1.0;
    }
    ux /= ul;
    uy /= ul;

    const outerHere = shapeField?.outer ? _sampleFieldClamped(shapeField.outer, W, H, p[0], p[1]) : 0.0;
    const coreHere = shapeField?.core ? _sampleFieldClamped(shapeField.core, W, H, p[0], p[1]) : 0.0;
    const contourHere = shapeField?.contour ? _sampleFieldClamped(shapeField.contour, W, H, p[0], p[1]) : 0.0;
    const out1 = shapeField?.support ? _sampleFieldClamped(shapeField.support, W, H, p[0] + ux * outwardPx1, p[1] + uy * outwardPx1) : 0.0;
    const out2 = shapeField?.support ? _sampleFieldClamped(shapeField.support, W, H, p[0] + ux * outwardPx2, p[1] + uy * outwardPx2) : 0.0;
    const outwardOpen = 1.0 - Math.max(out1, 0.72 * out2);

    const terminalRaw =
      0.40 * _smoothstep(0.50, 0.96, pr1) +
      0.32 * _smoothstep(0.12, 0.82, graphEnd) +
      0.18 * _smoothstep(0.34, 0.96, radial) +
      0.10 * _smoothstep(0.08, 0.70, outwardOpen);
    const terminalK = _clamp01(terminalRaw * (0.34 + 0.66 * longDominance) * (0.78 + 0.22 * outerHere) * (1.0 - 0.34 * coreHere));

    const cornerK = _clamp01((0.62 * bend + 0.24 * contourHere + 0.14 * graphEnd) * (0.40 + 0.60 * terminalK));
    const branchK = _clamp01((0.50 * Math.max(graphBranch, topoBranch) + 0.22 * contourHere + 0.12 * outwardOpen) * (0.26 + 0.74 * terminalK));
    const peripheralK = _clamp01(0.46 * outerHere + 0.24 * contourHere + 0.18 * radial + 0.12 * outwardOpen);
    const salientK = _clamp01(0.58 * terminalK + 0.16 * cornerK + 0.10 * branchK + 0.16 * peripheralK);

    terminal[i] = terminalK;
    corner[i] = cornerK;
    branch[i] = branchK;
    peripheral[i] = peripheralK;
    salient[i] = salientK;
  }

  return { terminal, corner, branch, peripheral, salient, graph, axes };
}


function _buildAFMHeightMergeLayer(W, H, coords, atoms, bonds, z_view, opts = {}) {
  const merge = new Float32Array(W * H);
  const mergeMask = new Float32Array(W * H);
  const peakDamp = new Float32Array(W * H);
  if (!bonds || !bonds.length) return { merge, mergeMask, peakDamp };

  const hide_front = !!opts.hide_front;
  const focal_z = Number.isFinite(opts.focal_z) ? opts.focal_z : 0.0;
  const atom_size_mul = Number.isFinite(opts.atom_size_mul) ? opts.atom_size_mul : 1.2;
  const atom_size_exp = Number.isFinite(opts.atom_size_exp) ? opts.atom_size_exp : 1.25;
  const scale = Number.isFinite(opts.scale) ? opts.scale : 1.0;

  const adj = _afmBuildBondAdjacency(atoms.length, bonds);
  const linearity = _afmAtomLinearityScores(coords, adj);

  for (let k = 0; k < bonds.length; k++) {
    const b = bonds[k];
    const i = b[0] | 0;
    const j = b[1] | 0;
    const a1 = atoms[i];
    const a2 = atoms[j];
    const p1 = coords[i];
    const p2 = coords[j];
    if (!a1 || !a2 || !p1 || !p2) continue;

    const z1 = z_view ? (z_view[i] ?? 0) : (a1.z ?? 0);
    const z2 = z_view ? (z_view[j] ?? 0) : (a2.z ?? 0);
    if (hide_front && (z1 > focal_z + 1e-9 || z2 > focal_z + 1e-9)) continue;

    const Z1 = a1.Z | 0;
    const Z2 = a2.Z | 0;
    const r1 = get_covalent_radius(Z1);
    const r2 = get_covalent_radius(Z2);
    let s1 = Math.max(1e-6, Math.pow(r1, atom_size_exp) * scale * atom_size_mul);
    let s2 = Math.max(1e-6, Math.pow(r2, atom_size_exp) * scale * atom_size_mul);

    const L = Math.hypot(p2[0] - p1[0], p2[1] - p1[1]);
    if (!(L > 1e-6)) continue;

    const sigmaLine = Math.max(0.8, Math.min(18.0, 0.62 * Math.min(s1, s2)));
    const bridgeStrengthBase = 0.18 + 0.14 * Math.max(linearity[i] || 0, linearity[j] || 0);
    const degreeBoost = 0.04 * Math.min(3, Math.max((adj[i]?.length || 0) - 1, 0) + Math.max((adj[j]?.length || 0) - 1, 0));
    const shortBondBoost = 0.08 * _smoothstep(2.2 * sigmaLine, 0.9 * sigmaLine, L);
    const bridgeStrength = Math.min(0.46, bridgeStrengthBase + degreeBoost + shortBondBoost);
    if (!(bridgeStrength > 1e-5)) continue;

    const maskStrength = Math.min(1.0, 0.55 + 0.45 * (bridgeStrength / 0.46));
    const dampStrength = Math.min(0.72, 0.26 + 0.46 * (bridgeStrength / 0.46));

    draw_gaussian_line({ a: merge, w: W, h: H, bg: 0 }, p1, p2, bridgeStrength, sigmaLine);
    draw_gaussian_line({ a: mergeMask, w: W, h: H, bg: 0 }, p1, p2, maskStrength, Math.max(0.9, sigmaLine * 1.05));
    draw_gaussian_line({ a: peakDamp, w: W, h: H, bg: 0 }, p1, p2, dampStrength, Math.max(0.8, sigmaLine * 0.9));
  }

  let mmax = 0.0;
  let maskMax = 0.0;
  let dampMax = 0.0;
  for (let i = 0; i < merge.length; i++) {
    if (merge[i] > mmax) mmax = merge[i];
    if (mergeMask[i] > maskMax) maskMax = mergeMask[i];
    if (peakDamp[i] > dampMax) dampMax = peakDamp[i];
  }
  if (mmax > 1e-12) {
    const inv = 1.0 / mmax;
    for (let i = 0; i < merge.length; i++) merge[i] *= inv;
  }
  if (maskMax > 1e-12) {
    const inv = 1.0 / maskMax;
    for (let i = 0; i < mergeMask.length; i++) mergeMask[i] *= inv;
  }
  if (dampMax > 1e-12) {
    const inv = 1.0 / dampMax;
    for (let i = 0; i < peakDamp.length; i++) peakDamp[i] *= inv;
  }

  return { merge, mergeMask, peakDamp };
}

function _applyAFMLocalHeightMerge(height, W, H, mergeData) {
  if (!mergeData) return { mergePresence: new Float32Array(W * H), peakSoft: new Float32Array(W * H) };
  const merge = mergeData.merge || new Float32Array(W * H);
  const mergeMask = mergeData.mergeMask || new Float32Array(W * H);
  const peakDamp = mergeData.peakDamp || new Float32Array(W * H);

  const soft = new Float32Array(height);
  const scratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: soft, w: W, h: H, bg: 0 }, 1.2, scratch);

  const mergePresence = new Float32Array(W * H);
  const peakSoft = new Float32Array(W * H);

  for (let p = 0; p < height.length; p++) {
    const m = merge[p];
    const mm = mergeMask[p];
    const pd = peakDamp[p];
    const h = height[p];
    const s = soft[p];

    const saddleFill = Math.max(h, s * (0.62 + 0.28 * m));
    const peakBlend = Math.min(0.34, 0.18 * mm + 0.16 * pd);
    const merged = h * (1.0 - peakBlend) + saddleFill * peakBlend;
    const localMask = _clamp01(Math.max(mm * 0.92, m * 0.82) * _smoothstep(0.05, 0.32, Math.max(h, s)));

    height[p] = merged;
    mergePresence[p] = localMask;
    peakSoft[p] = _clamp01(pd * _smoothstep(0.12, 0.70, h));
  }

  let hmax = 0.0;
  for (let p = 0; p < height.length; p++) if (height[p] > hmax) hmax = height[p];
  if (hmax > 1e-12) {
    const inv = 1.0 / hmax;
    for (let p = 0; p < height.length; p++) height[p] *= inv;
  }

  return { mergePresence, peakSoft };
}


function _buildAFMBondBlurCompLayers(bondLayer, W, H, blurSigma, waveWidthPx, ctx = {}) {
  const blurK = _clamp01((Math.max(0, blurSigma) - 1.0) / 3.4);
  const out = {
    core: bondLayer.a,
    local: null,
    wide: null,
    gate: null,
    open: null,
    contour: null,
    blurK,
  };
  if (!(blurK > 1e-6)) return out;

  const sigmaNear = Math.max(
    0.6,
    Math.min(2.7, 0.46 * Math.max(blurSigma, 0) + 0.08 * Math.max(waveWidthPx || 0, 0)),
  );
  const sigmaWide = Math.max(sigmaNear + 0.35, Math.min(5.2, sigmaNear * 1.95));

  const nearSeed = new Float32Array(W * H);
  const wideSeed = new Float32Array(W * H);
  const gate = new Float32Array(W * H);
  const openArr = new Float32Array(W * H);
  const contourArr = new Float32Array(W * H);

  for (let p = 0; p < bondLayer.a.length; p++) {
    const b = bondLayer.a[p];
    if (!(b > 1e-8)) continue;

    const h = ctx.height ? ctx.height[p] : 0.0;
    const merge = ctx.mergePresence ? ctx.mergePresence[p] : 0.0;
    const mol = ctx.molMask ? ctx.molMask[p] : 0.0;
    const bondMask = ctx.bondMask ? ctx.bondMask[p] : 0.0;
    const feature = ctx.bondFeature ? ctx.bondFeature[p] : 0.0;
    const openK = ctx.openField ? ctx.openField[p] : 0.0;
    const contourK = ctx.contourField ? ctx.contourField[p] : 0.0;
    const endpointK = ctx.endpointField ? ctx.endpointField[p] : 0.0;

    const relief = _smoothstep(0.04, 0.46, h);
    const support = _clamp01(Math.max(mol, 0.78 * relief, 0.72 * merge, 0.62 * bondMask));
    const edgeZone = _clamp01(Math.max(contourK, 0.70 * merge, 0.56 * bondMask));
    const openZone = _clamp01(Math.max(openK, endpointK, 0.78 * feature));
    const dense = _smoothstep(0.48, 0.92, Math.max(h, mol)) * (1.0 - 0.34 * merge);

    const gateP = _clamp01(
      0.08 +
      0.80 * support +
      0.24 * edgeZone +
      0.16 * openZone -
      0.10 * (1.0 - support) * (1.0 - edgeZone),
    );

    gate[p] = gateP;
    openArr[p] = openZone;
    contourArr[p] = edgeZone;

    nearSeed[p] =
      b *
      gateP *
      (0.82 + 0.18 * feature) *
      (1.0 - 0.16 * dense + 0.08 * edgeZone);

    wideSeed[p] =
      b *
      gateP *
      (0.10 + 0.52 * openZone + 0.20 * edgeZone + 0.18 * feature) *
      (1.0 - 0.56 * dense);
  }

  const scratchNear = new Float32Array(W * H);
  const local = new Float32Array(nearSeed);
  gaussianBlurFloat({ a: local, w: W, h: H, bg: 0 }, sigmaNear, scratchNear);

  const scratchWide = new Float32Array(W * H);
  const wide = new Float32Array(wideSeed);
  gaussianBlurFloat({ a: wide, w: W, h: H, bg: 0 }, sigmaWide, scratchWide);

  for (let p = 0; p < local.length; p++) {
    const gateP = gate[p];
    if (!(gateP > 1e-8)) {
      local[p] = 0.0;
      wide[p] = 0.0;
      continue;
    }

    local[p] *= gateP;
    const wideGate = _clamp01(
      0.04 +
      0.56 * gateP +
      0.26 * openArr[p] +
      0.16 * contourArr[p],
    );
    wide[p] *= wideGate * (1.0 - 0.42 * (1.0 - gateP));
  }

  out.local = local;
  out.wide = wide;
  out.gate = gate;
  out.open = openArr;
  out.contour = contourArr;
  return out;
}

function _composeAFMBondValue(p, bondComp, ctx = {}) {
  const core = bondComp?.core ? bondComp.core[p] : 0.0;
  const local = bondComp?.local ? bondComp.local[p] : 0.0;
  const wide = bondComp?.wide ? bondComp.wide[p] : 0.0;
  if (!(core > 1e-8 || local > 1e-8 || wide > 1e-8)) return 0.0;

  const blurK = bondComp?.blurK || 0.0;
  const h = ctx.height ? ctx.height[p] : 0.0;
  const merge = ctx.mergePresence ? ctx.mergePresence[p] : 0.0;
  const mol = ctx.molMask ? ctx.molMask[p] : 0.0;
  const feature = ctx.bondFeature ? ctx.bondFeature[p] : 0.0;
  const support = bondComp?.gate ? bondComp.gate[p] : 0.0;
  const openK = bondComp?.open ? bondComp.open[p] : 0.0;
  const contourK = bondComp?.contour ? bondComp.contour[p] : 0.0;

  const dense = _smoothstep(0.50, 0.93, Math.max(h, mol)) * (1.0 - 0.24 * merge);
  const coreGain = 1.0 + 0.05 * blurK + 0.08 * feature;
  const localGain =
    blurK *
    (0.08 + 0.22 * support + 0.20 * contourK + 0.22 * openK) *
    (1.0 - 0.22 * dense);
  const wideGain =
    blurK *
    blurK *
    (0.008 + 0.070 * openK + 0.050 * contourK + 0.022 * feature) *
    (1.0 - 0.34 * dense);

  return core * coreGain + local * localGain + wide * wideGain;
}

function _buildAFMBondTopologyField(W, H, coords, atoms, bonds, z_view, opts = {}) {
  const feature = new Float32Array(W * H);
  const tipField = new Float32Array(W * H);
  const bendField = new Float32Array(W * H);
  const endpointField = new Float32Array(W * H);
  const cornerField = new Float32Array(W * H);
  const branchField = new Float32Array(W * H);
  const terminalField = new Float32Array(W * H);
  const terminalSpillField = new Float32Array(W * H);
  const openField = new Float32Array(W * H);
  const contourField = new Float32Array(W * H);
  const mask = new Float32Array(W * H);

  if (!bonds || !bonds.length) {
    return {
      feature,
      tipField,
      bendField,
      endpointField,
      cornerField,
      branchField,
      terminalField,
      terminalSpillField,
      openField,
      contourField,
      mask,
      topology: _afmAtomTopologyScores(coords, _afmBuildBondAdjacency(atoms.length, bonds || [])),
    };
  }

  const hide_front = !!opts.hide_front;
  const focal_z = Number.isFinite(opts.focal_z) ? opts.focal_z : 0.0;
  const atom_size_mul = Number.isFinite(opts.atom_size_mul) ? opts.atom_size_mul : 1.2;
  const atom_size_exp = Number.isFinite(opts.atom_size_exp) ? opts.atom_size_exp : 1.25;
  const scale = Number.isFinite(opts.scale) ? opts.scale : 1.0;

  const ctx = {
    height: opts.height || null,
    mergePresence: opts.mergePresence || null,
    molMask: opts.molMask || null,
    shapeSupport: opts.shapeSupport || null,
    shapeContour: opts.shapeContour || null,
    shapeOuter: opts.shapeOuter || null,
    shapeCore: opts.shapeCore || null,
  };

  const adj = _afmBuildBondAdjacency(atoms.length, bonds);
  const topo = _afmAtomTopologyScores(coords, adj);
  const coarse = _buildAFMCoarseShapeField(coords, atoms, adj, topo, {
    support: ctx.shapeSupport,
    contour: ctx.shapeContour,
    outer: ctx.shapeOuter,
    core: ctx.shapeCore,
  }, W, H, Number.isFinite(opts.angstroms_per_pixel) ? opts.angstroms_per_pixel : 0.1);

  for (let k = 0; k < bonds.length; k++) {
    const b = bonds[k];
    const i = b[0] | 0;
    const j = b[1] | 0;
    const a1 = atoms[i];
    const a2 = atoms[j];
    const p1 = coords[i];
    const p2 = coords[j];
    if (!a1 || !a2 || !p1 || !p2) continue;

    const z1 = z_view ? (z_view[i] ?? 0) : (a1.z ?? 0);
    const z2 = z_view ? (z_view[j] ?? 0) : (a2.z ?? 0);
    if (hide_front && (z1 > focal_z + 1e-9 || z2 > focal_z + 1e-9)) continue;

    const Z1 = a1.Z | 0;
    const Z2 = a2.Z | 0;
    const r1 = get_covalent_radius(Z1);
    const r2 = get_covalent_radius(Z2);
    const s1 = Math.max(1e-6, Math.pow(r1, atom_size_exp) * scale * atom_size_mul);
    const s2 = Math.max(1e-6, Math.pow(r2, atom_size_exp) * scale * atom_size_mul);
    const sigmaLine = Math.max(0.75, Math.min(18.0, 0.58 * Math.min(s1, s2)));

    const tipI = topo.tip[i] || 0.0;
    const tipJ = topo.tip[j] || 0.0;
    const bendI = topo.bend[i] || 0.0;
    const bendJ = topo.bend[j] || 0.0;
    const branchI = topo.branch[i] || 0.0;
    const branchJ = topo.branch[j] || 0.0;
    const expI = topo.exposure[i] || 0.0;
    const expJ = topo.exposure[j] || 0.0;
    const termI = coarse.terminal[i] || 0.0;
    const termJ = coarse.terminal[j] || 0.0;
    const coarseCornerI = coarse.corner[i] || 0.0;
    const coarseCornerJ = coarse.corner[j] || 0.0;
    const coarseBranchI = coarse.branch[i] || 0.0;
    const coarseBranchJ = coarse.branch[j] || 0.0;
    const coarsePeripheralI = coarse.peripheral[i] || 0.0;
    const coarsePeripheralJ = coarse.peripheral[j] || 0.0;
    const coarseSalientI = coarse.salient[i] || 0.0;
    const coarseSalientJ = coarse.salient[j] || 0.0;

    const endpointTopo = Math.max(tipI, tipJ, 0.86 * Math.max(bendI, bendJ));
    const lineFeature = _clamp01(
      0.18 * Math.max(expI, expJ) +
      0.22 * Math.max(branchI, branchJ) +
      0.60 * endpointTopo,
    );
    let edgeOpen = 0.0;
    let endOpen1 = 0.0;
    let endOpen2 = 0.0;
    let midDense = 0.0;
    let openTrust = 1.0;
    let shapeContourAvg = 0.0;
    let shapeOuterAvg = 0.0;
    let shapeCoreAvg = 0.0;

    const dx = p2[0] - p1[0];
    const dy = p2[1] - p1[1];
    const L = Math.hypot(dx, dy);
    if (L > 1e-6 && (ctx.height || ctx.molMask || ctx.mergePresence)) {
      const tx = dx / L;
      const ty = dy / L;
      const nx = -ty;
      const ny = tx;
      const sideDist = Math.max(1.2, sigmaLine * 1.18);
      const sideDist2 = Math.max(sideDist + 0.6, sigmaLine * 1.75);
      const endDist = Math.max(1.3, sigmaLine * 1.20);
      const edgeMargin = Math.max(2.5, sideDist2 + sigmaLine * 0.9);

      const ts = [0.22, 0.50, 0.78];
      let openSum = 0.0;
      let denseSum = 0.0;
      let trustSum = 0.0;
      let contourSum = 0.0;
      let outerSum = 0.0;
      let coreSum = 0.0;
      let cnt = 0;

      for (let s = 0; s < ts.length; s++) {
        const t = ts[s];
        const x = p1[0] + dx * t;
        const y = p1[1] + dy * t;

        const c = _sampleAFMStructureSupport(ctx, W, H, x, y);
        const left1 = _sampleAFMStructureSupportSafe(ctx, W, H, x + nx * sideDist, y + ny * sideDist, edgeMargin, c);
        const right1 = _sampleAFMStructureSupportSafe(ctx, W, H, x - nx * sideDist, y - ny * sideDist, edgeMargin, c);
        const left2 = _sampleAFMStructureSupportSafe(ctx, W, H, x + nx * sideDist2, y + ny * sideDist2, edgeMargin, c);
        const right2 = _sampleAFMStructureSupportSafe(ctx, W, H, x - nx * sideDist2, y - ny * sideDist2, edgeMargin, c);

        const left = Math.max(left1.value, 0.72 * left2.value);
        const right = Math.max(right1.value, 0.72 * right2.value);
        const sideTrust = Math.min(left1.conf, right1.conf, left2.conf, right2.conf);

        const sideGap = 1.0 - Math.min(left, right);
        const sideAsym = Math.abs(left - right);
        const openHere = _clamp01(0.72 * sideGap + 0.28 * sideAsym) * sideTrust;

        openSum += openHere;
        denseSum += Math.max(c, 0.65 * (left + right) * 0.5);
        trustSum += sideTrust;
        if (ctx.shapeContour) contourSum += _sampleFieldBilinear(ctx.shapeContour, W, H, x, y);
        if (ctx.shapeOuter) outerSum += _sampleFieldBilinear(ctx.shapeOuter, W, H, x, y);
        if (ctx.shapeCore) coreSum += _sampleFieldBilinear(ctx.shapeCore, W, H, x, y);
        cnt++;
      }

      if (cnt > 0) {
        const trustAvg = trustSum / cnt;
        edgeOpen = openSum / cnt;
        midDense = denseSum / cnt;
        openTrust = trustAvg;
        shapeContourAvg = contourSum / cnt;
        shapeOuterAvg = outerSum / cnt;
        shapeCoreAvg = coreSum / cnt;
      }

      const endMargin = Math.max(2.6, endDist * 1.8 + sigmaLine * 0.8);
      const p1Back1 = _sampleAFMStructureSupportSafe(ctx, W, H, p1[0] - tx * endDist, p1[1] - ty * endDist, endMargin, 1.0);
      const p1Back2 = _sampleAFMStructureSupportSafe(ctx, W, H, p1[0] - tx * (endDist * 1.8), p1[1] - ty * (endDist * 1.8), endMargin, 1.0);
      const p2Front1 = _sampleAFMStructureSupportSafe(ctx, W, H, p2[0] + tx * endDist, p2[1] + ty * endDist, endMargin, 1.0);
      const p2Front2 = _sampleAFMStructureSupportSafe(ctx, W, H, p2[0] + tx * (endDist * 1.8), p2[1] + ty * (endDist * 1.8), endMargin, 1.0);

      const endTrust1 = Math.min(p1Back1.conf, p1Back2.conf);
      const endTrust2 = Math.min(p2Front1.conf, p2Front2.conf);

      endOpen1 = _clamp01(
        1.0 - Math.max(p1Back1.value, 0.72 * p1Back2.value),
      ) * endTrust1;
      endOpen2 = _clamp01(
        1.0 - Math.max(p2Front1.value, 0.72 * p2Front2.value),
      ) * endTrust2;

      openTrust = Math.min(1.0, Math.max(openTrust, 0.5 * (endTrust1 + endTrust2)));
    }


    const coarseEnd = Math.max(endOpen1, endOpen2);
    const terminalSeed = Math.max(termI, termJ, 0.82 * Math.max(coarseSalientI, coarseSalientJ));
    const midx = 0.5 * (p1[0] + p2[0]);
    const midy = 0.5 * (p1[1] + p2[1]);
    const axisDx = midx - coarse.axes.cx;
    const axisDy = midy - coarse.axes.cy;
    const longPos = Math.abs(axisDx * coarse.axes.ax1x + axisDy * coarse.axes.ax1y) / Math.max(1e-6, coarse.axes.maxProj1);
    const crossPos = Math.abs(axisDx * coarse.axes.ax2x + axisDy * coarse.axes.ax2y) / Math.max(1e-6, coarse.axes.maxProj2);
    const longEndBand = _smoothstep(0.56, 0.90, longPos);
    const crossKeep = 1.0 - 0.22 * _smoothstep(0.82, 1.06, crossPos);
    const terminalRegion = _clamp01(
      longEndBand *
      crossKeep *
      (0.34 + 0.28 * Math.max(shapeOuterAvg, Math.max(coarsePeripheralI, coarsePeripheralJ)) + 0.18 * edgeOpen + 0.20 * terminalSeed) *
      (1.0 - 0.14 * shapeCoreAvg)
    );
    const terminalGate = _clamp01(
      0.28 * terminalSeed +
      0.18 * coarseEnd +
      0.10 * shapeOuterAvg +
      0.08 * Math.max(coarsePeripheralI, coarsePeripheralJ) +
      0.04 * edgeOpen +
      0.56 * terminalRegion -
      0.10 * shapeCoreAvg,
    );
    const endpointAmp = _clamp01(
      (0.10 * Math.max(tipI, tipJ) + 0.62 * coarseEnd + 0.30 * terminalSeed + 0.30 * terminalRegion + 0.08 * shapeOuterAvg) *
      openTrust *
      (0.46 + 0.54 * terminalGate),
    );
    const cornerAmp = _clamp01(
      (0.16 * Math.max(bendI, bendJ) + 0.34 * Math.max(coarseCornerI, coarseCornerJ) + 0.18 * shapeContourAvg + 0.22 * terminalRegion) *
      (0.34 + 0.66 * edgeOpen) *
      (0.54 + 0.46 * terminalGate) *
      (0.72 + 0.28 * openTrust),
    );
    const branchAmp = _clamp01(
      (0.16 * Math.max(branchI, branchJ) + 0.32 * Math.max(coarseBranchI, coarseBranchJ) + 0.12 * shapeContourAvg + 0.10 * shapeOuterAvg + 0.24 * terminalRegion) *
      (0.30 + 0.70 * edgeOpen) *
      (0.50 + 0.50 * terminalGate) *
      (0.78 + 0.22 * openTrust),
    );
    const openAmp = _clamp01(
      (
        0.16 * edgeOpen +
        0.18 * coarseEnd +
        0.20 * terminalSeed +
        0.28 * terminalRegion +
        0.10 * shapeOuterAvg +
        0.08 * Math.max(expI, expJ)
      ) * openTrust * (0.42 + 0.58 * terminalGate),
    );
    const contourAmp = _clamp01(
      (
        0.18 * shapeContourAvg +
        0.14 * terminalSeed +
        0.26 * terminalRegion +
        0.16 * openAmp +
        0.14 * cornerAmp +
        0.06 * branchAmp +
        0.06 * Math.max(0.0, shapeOuterAvg - 0.6 * shapeCoreAvg)
      ) * (0.42 + 0.58 * terminalGate),
    );

    const denseDamp = 1.0 - 0.34 * _smoothstep(0.48, 0.88, midDense) * (1.0 - 0.78 * openAmp);
    const selective = _clamp01(
      (
        0.03 * lineFeature +
        0.30 * endpointAmp +
        0.18 * cornerAmp +
        0.08 * branchAmp +
        0.14 * openAmp +
        0.12 * contourAmp +
        0.15 * terminalRegion
      ) * denseDamp * (0.52 + 0.48 * terminalGate),
    );

    draw_gaussian_line(
      { a: feature, w: W, h: H, bg: 0 },
      p1,
      p2,
      0.12 + 0.62 * selective,
      Math.max(0.8, sigmaLine * 0.96),
    );
    draw_gaussian_line(
      { a: terminalField, w: W, h: H, bg: 0 },
      p1,
      p2,
      0.10 + 0.42 * terminalGate + 0.46 * terminalRegion,
      Math.max(0.95, sigmaLine * 1.18),
    );
    draw_gaussian_line(
      { a: openField, w: W, h: H, bg: 0 },
      p1,
      p2,
      0.10 + 0.52 * openAmp,
      Math.max(0.8, sigmaLine * 1.02),
    );
    draw_gaussian_line(
      { a: contourField, w: W, h: H, bg: 0 },
      p1,
      p2,
      0.10 + 0.52 * contourAmp,
      Math.max(0.9, sigmaLine * 1.18),
    );
    draw_gaussian_line(
      { a: mask, w: W, h: H, bg: 0 },
      p1,
      p2,
      1.0,
      Math.max(0.9, sigmaLine * 1.18),
    );

    const endpointSigma = Math.max(0.95, sigmaLine * 1.14);
    if (tipI > 1e-6) _drawAFMGaussianSpot(tipField, W, H, p1[0], p1[1], 0.46 + 0.54 * tipI, endpointSigma);
    if (tipJ > 1e-6) _drawAFMGaussianSpot(tipField, W, H, p2[0], p2[1], 0.46 + 0.54 * tipJ, endpointSigma);
    if (bendI > 1e-6) _drawAFMGaussianSpot(bendField, W, H, p1[0], p1[1], 0.34 + 0.66 * bendI, endpointSigma * 0.95);
    if (bendJ > 1e-6) _drawAFMGaussianSpot(bendField, W, H, p2[0], p2[1], 0.34 + 0.66 * bendJ, endpointSigma * 0.95);

    if (endpointAmp > 1e-6) {
      _drawAFMGaussianSpot(endpointField, W, H, p1[0], p1[1], 0.16 + 0.84 * Math.max(endpointAmp, tipI), endpointSigma * 1.02);
      _drawAFMGaussianSpot(endpointField, W, H, p2[0], p2[1], 0.16 + 0.84 * Math.max(endpointAmp, tipJ), endpointSigma * 1.02);
    }
    if (cornerAmp > 1e-6) {
      _drawAFMGaussianSpot(cornerField, W, H, p1[0], p1[1], 0.12 + 0.88 * Math.max(cornerAmp, bendI), endpointSigma * 0.98);
      _drawAFMGaussianSpot(cornerField, W, H, p2[0], p2[1], 0.12 + 0.88 * Math.max(cornerAmp, bendJ), endpointSigma * 0.98);
    }
    if (branchAmp > 1e-6) {
      _drawAFMGaussianSpot(branchField, W, H, p1[0], p1[1], 0.10 + 0.90 * Math.max(branchAmp, branchI), endpointSigma * 1.02);
      _drawAFMGaussianSpot(branchField, W, H, p2[0], p2[1], 0.10 + 0.90 * Math.max(branchAmp, branchJ), endpointSigma * 1.02);
    }
  }

  const lineN = _normalizePositive(feature, 0.96);
  const tipN = _normalizePositive(tipField, 0.92);
  const bendN = _normalizePositive(bendField, 0.95);
  const endpointN = _normalizePositive(endpointField, 0.95);
  const cornerN = _normalizePositive(cornerField, 0.95);
  const branchN = _normalizePositive(branchField, 0.95);
  const terminalN = _normalizePositive(terminalField, 0.94);
  const terminalWide = new Float32Array(terminalN);
  const terminalWideScratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: terminalWide, w: W, h: H, bg: 0 }, 1.35, terminalWideScratch);
  const terminalWideN = _normalizePositive(terminalWide, 0.96);
  const terminalSpill = new Float32Array(terminalN);
  const terminalSpillScratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: terminalSpill, w: W, h: H, bg: 0 }, 2.25, terminalSpillScratch);
  const terminalSpillN = _normalizePositive(terminalSpill, 0.98);
  const openN = _normalizePositive(openField, 0.96);
  const contourN = _normalizePositive(contourField, 0.98);
  const maskN = _normalizePositive(mask, 1.0);

  for (let p = 0; p < feature.length; p++) {
    const endpointAct = _smoothstep(0.16, 0.62, Math.max(endpointN[p], tipN[p]));
    const cornerAct = _smoothstep(0.18, 0.64, Math.max(cornerN[p], 0.84 * bendN[p]));
    const terminalCoreAct = _smoothstep(0.14, 0.56, terminalN[p]);
    const terminalRegionAct = _smoothstep(0.10, 0.48, terminalWideN[p]);
    const terminalSpillAct = _smoothstep(0.08, 0.42, terminalSpillN[p]);
    const terminalAct = _clamp01(0.42 * terminalCoreAct + 0.58 * terminalRegionAct);
    const terminalEnvAct = _clamp01(Math.max(terminalAct, 0.68 * terminalSpillAct));
    const openAct = _smoothstep(0.18, 0.66, openN[p]) * (0.66 + 0.34 * terminalEnvAct);
    const contourAct = _smoothstep(0.20, 0.68, contourN[p]) * (0.64 + 0.36 * terminalEnvAct);
    const branchAct = _smoothstep(0.20, 0.68, branchN[p]) * (0.70 + 0.30 * terminalEnvAct);
    const lineAct = _smoothstep(0.22, 0.72, lineN[p]);

    feature[p] = _clamp01(
      0.06 * lineAct +
      0.26 * endpointAct +
      0.16 * cornerAct +
      0.08 * branchAct +
      0.16 * openAct +
      0.18 * contourAct +
      0.08 * terminalAct +
      0.08 * terminalEnvAct,
    );
    tipField[p] = tipN[p];
    bendField[p] = bendN[p];
    endpointField[p] = endpointAct;
    cornerField[p] = cornerAct;
    branchField[p] = branchAct;
    terminalField[p] = terminalAct;
    terminalSpillField[p] = terminalSpillAct;
    openField[p] = openAct;
    contourField[p] = contourAct;
    mask[p] = maskN[p];
  }

  return {
    feature,
    tipField,
    bendField,
    endpointField,
    cornerField,
    branchField,
    terminalField,
    terminalSpillField,
    openField,
    contourField,
    mask,
    topology: topo,
  };
}

function _buildMoleculeRadialField(height, W, H, coords, atoms, angstroms_per_pixel, opts = {}) {
  const sigma = Math.max(3.0, Math.min(24.0, 0.035 * Math.min(W, H) + 4.0));
  const broad = new Float32Array(height);
  const scratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: broad, w: W, h: H, bg: 0 }, sigma, scratch);
  const broadN = _normalizePositive(broad, 0.88);

  let sw = 0.0;
  let sx = 0.0;
  let sy = 0.0;
  for (let i = 0; i < atoms.length; i++) {
    const p = coords[i];
    const a = atoms[i];
    if (!p || !a) continue;
    const Z = (a.Z | 0) || 6;
    const w = Math.max(0.25, get_covalent_radius(Z));
    sw += w;
    sx += p[0] * w;
    sy += p[1] * w;
  }

  if (!(sw > 1e-8)) {
    for (let y = 0; y < H; y++) {
      const row = y * W;
      for (let x = 0; x < W; x++) {
        const w = Math.pow(Math.max(0.0, broadN[row + x]), 1.35);
        if (w > 1e-8) {
          sw += w;
          sx += x * w;
          sy += y * w;
        }
      }
    }
  }

  const cx = sw > 1e-8 ? sx / sw : 0.5 * (W - 1);
  const cy = sw > 1e-8 ? sy / sw : 0.5 * (H - 1);

  const radial = new Float32Array(W * H);
  const radialBond = new Float32Array(W * H);
  const radialAtom = new Float32Array(W * H);
  const radialBlur = new Float32Array(W * H);
  const molMask = new Float32Array(W * H);
  const distAField = new Float32Array(W * H);

  for (let y = 0; y < H; y++) {
    const row = y * W;
    for (let x = 0; x < W; x++) {
      const p = row + x;
      const dx = (x - cx) * angstroms_per_pixel;
      const dy = (y - cy) * angstroms_per_pixel;
      const distA = Math.hypot(dx, dy);
      const mol = _smoothstep(0.02, 0.16, broadN[p]);
      const prof = _radialProfileValue(distA, opts.profile);
      const k = prof * mol;
      distAField[p] = distA;
      molMask[p] = mol;
      radial[p] = k;
      radialBond[p] = k;
      radialAtom[p] = Math.pow(k, 1.22);
      radialBlur[p] = Math.pow(k, 0.95) * mol;
    }
  }

  return { cx, cy, broadN, radial, radialBond, radialAtom, radialBlur, molMask, distAField };
}

function _buildAFMGlobalShapeField(height, W, H, mergeData, mergeState, angstroms_per_pixel) {
  const merge = mergeData?.merge || new Float32Array(W * H);
  const mergeMask = mergeData?.mergeMask || new Float32Array(W * H);
  const mergePresence = mergeState?.mergePresence || new Float32Array(W * H);

  const base = new Float32Array(W * H);
  for (let p = 0; p < base.length; p++) {
    base[p] = Math.max(
      0.96 * mergeMask[p],
      0.82 * merge[p],
      0.68 * mergePresence[p],
      0.52 * height[p],
    );
  }

  const sigmaNear = Math.max(2.0, Math.min(20.0, 1.65 / Math.max(angstroms_per_pixel, 1e-6)));
  const sigmaWide = Math.max(sigmaNear + 1.2, Math.min(34.0, 3.05 / Math.max(angstroms_per_pixel, 1e-6)));

  const near = new Float32Array(base);
  const nearScratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: near, w: W, h: H, bg: 0 }, sigmaNear, nearScratch);

  const wide = new Float32Array(base);
  const wideScratch = new Float32Array(W * H);
  gaussianBlurFloat({ a: wide, w: W, h: H, bg: 0 }, sigmaWide, wideScratch);

  const nearN = _normalizePositive(near, 0.96);
  const wideN = _normalizePositive(wide, 0.98);

  const support = new Float32Array(W * H);
  const core = new Float32Array(W * H);
  const contour = new Float32Array(W * H);
  const outer = new Float32Array(W * H);

  for (let p = 0; p < support.length; p++) {
    const n = nearN[p];
    const w = wideN[p];
    const marker = Math.max(n, 0.86 * w);
    const coreK = _smoothstep(0.22, 0.54, n);
    const supportK = _smoothstep(0.10, 0.32, marker);
    const contourRaw = Math.max(0.0, n - 0.74 * w);
    const contourK = _smoothstep(0.05, 0.22, contourRaw * 2.8);
    const outerRaw = Math.max(0.0, supportK - 0.72 * coreK + 0.34 * contourK);
    const outerK = _smoothstep(0.06, 0.28, outerRaw);

    support[p] = supportK;
    core[p] = coreK;
    contour[p] = contourK;
    outer[p] = outerK;
  }

  return { base, near, wide, nearN, wideN, support, core, contour, outer };
}

function _applyAFMContourShaping(height, W, H, mergeData, mergeState, bondTopology, radialField, globalShapeField) {
  if (!bondTopology || !bondTopology.contourField) return;

  const contour = bondTopology.contourField;
  const openField = bondTopology.openField || new Float32Array(W * H);
  const endpointField = bondTopology.endpointField || new Float32Array(W * H);
  const feature = bondTopology.feature || new Float32Array(W * H);
  const mergePresence = mergeState?.mergePresence || new Float32Array(W * H);
  const peakSoft = mergeState?.peakSoft || new Float32Array(W * H);
  const mergeBase = mergeData?.merge || new Float32Array(W * H);
  const mergeMask = mergeData?.mergeMask || new Float32Array(W * H);
  const molMask = radialField?.molMask || new Float32Array(W * H);
  const shapeSupport = globalShapeField?.support || new Float32Array(W * H);
  const shapeCore = globalShapeField?.core || new Float32Array(W * H);
  const shapeContour = globalShapeField?.contour || new Float32Array(W * H);
  const shapeOuter = globalShapeField?.outer || new Float32Array(W * H);

  const softNear = new Float32Array(height);
  const scratchNear = new Float32Array(W * H);
  gaussianBlurFloat({ a: softNear, w: W, h: H, bg: 0 }, 0.95, scratchNear);

  const softWide = new Float32Array(height);
  const scratchWide = new Float32Array(W * H);
  gaussianBlurFloat({ a: softWide, w: W, h: H, bg: 0 }, 2.05, scratchWide);

  for (let p = 0; p < height.length; p++) {
    const contourK = contour[p];
    const openK = openField[p];
    const endpointK = endpointField[p];
    const featureK = feature[p];
    const mergeK = mergePresence[p];
    const peakK = peakSoft[p];
    const mergeLine = mergeBase[p];
    const mergeShape = mergeMask[p];
    const mol = molMask[p];
    const markerSupport = shapeSupport[p];
    const markerCore = shapeCore[p];
    const markerContour = shapeContour[p];
    const markerOuter = shapeOuter[p];

    const structK = Math.max(contourK, markerContour, 0.78 * mergeShape, 0.62 * mergeLine);
    if (!(structK > 1e-6)) continue;

    const h = height[p];
    const sNear = softNear[p];
    const sWide = softWide[p];
    const support = Math.max(markerSupport, mol, 0.76 * h, 0.72 * mergeK, 0.66 * mergeShape);

    const ridgeK = Math.max(structK, markerContour, 0.74 * mergeShape);
    const bridgeLift =
      0.12 * ridgeK +
      0.08 * mergeShape +
      0.06 * markerOuter +
      0.07 * openK +
      0.04 * endpointK +
      0.04 * featureK;
    const bridgeTarget = Math.max(
      h,
      sNear * (0.90 + 0.08 * support) +
        0.18 * ridgeK +
        0.10 * mergeK +
        0.10 * markerSupport +
        bridgeLift,
    );
    const bridgeMix = _clamp01(
      0.10 + 0.28 * ridgeK + 0.10 * openK + 0.08 * markerOuter + 0.06 * endpointK,
    );
    let shaped = h * (1.0 - bridgeMix) + bridgeTarget * bridgeMix;

    const capDelta = Math.max(0.0, shaped - sWide);
    const capFlatten =
      _clamp01(
        (0.10 + 0.28 * Math.max(contourK, markerContour) + 0.16 * Math.max(openK, markerOuter) + 0.06 * featureK) *
          (1.0 - 0.22 * mergeShape) *
          (1.0 - 0.32 * peakK),
      ) * _smoothstep(0.015, 0.12, capDelta);

    if (capFlatten > 1e-6) {
      const flattenTarget = Math.max(
        sNear * (0.98 + 0.02 * support),
        sWide * (0.98 + 0.04 * ridgeK),
        shaped * (0.92 + 0.04 * mergeK),
      );
      shaped = shaped * (1.0 - capFlatten) + flattenTarget * capFlatten;
    }

    const outerTighten =
      _clamp01(0.06 + 0.14 * Math.max(contourK, markerContour) + 0.10 * Math.max(openK, markerOuter)) *
      _smoothstep(0.02, 0.14, shaped - sWide) *
      (1.0 - 0.38 * mergeShape);
    if (outerTighten > 1e-6) {
      shaped = shaped * (1.0 - outerTighten) + Math.max(sWide, sNear * 0.98) * outerTighten;
    }

    height[p] = shaped;
  }

  let hmax = 0.0;
  for (let p = 0; p < height.length; p++) {
    if (height[p] > hmax) hmax = height[p];
  }
  if (hmax > 1e-12) {
    const inv = 1.0 / hmax;
    for (let p = 0; p < height.length; p++) height[p] *= inv;
  }
}


export function render_afm_like(atoms, opts = {}) {
  const ov = opts?.element_overrides || null;
  const bonds = "bonds" in opts ? opts.bonds : null;
  const img_size = Array.isArray(opts.img_size) ? opts.img_size : [400, 400];
  const angstroms_per_pixel = Number.isFinite(opts.angstroms_per_pixel)
    ? opts.angstroms_per_pixel
    : 0.1;
  const blur_sigma = Number.isFinite(opts.blur_sigma) ? opts.blur_sigma : 2.0;
  const background_gray = Number.isFinite(opts.background_gray)
    ? opts.background_gray
    : 127;
  const invert = !!opts.invert;
  const noise_stddev = Number.isFinite(opts.noise_stddev)
    ? opts.noise_stddev
    : 0.0;
  const contrast = Number.isFinite(opts.contrast) ? opts.contrast : 1.0;
  const compose_mode = opts.compose_mode || "sum"; // kept for API compat; not used here
  void compose_mode;

  const draw_bonds_flag = opts.draw_bonds_flag !== false;
  const camera = opts.camera || null;
  const bond_wave_width_px = Number.isFinite(opts.bond_wave_width_px)
    ? opts.bond_wave_width_px
    : 6.0;
  const bond_wave_amplitude = Number.isFinite(opts.bond_wave_amplitude)
    ? opts.bond_wave_amplitude
    : 0.4;

  const low_clip = opts.low_clip != null ? opts.low_clip : null;
  const high_clip = opts.high_clip != null ? opts.high_clip : null;

  const focal_z = Number.isFinite(opts.focal_z) ? opts.focal_z : 0.0;
  const hide_front = !!opts.hide_front;

  const show_scale_bar = !!opts.show_scale_bar;
  const scale_bar_corner = opts.scale_bar_corner || "bl";
  const scale_bar_margin_px = Number.isFinite(opts.scale_bar_margin_px)
    ? opts.scale_bar_margin_px
    : 12;

  // atom size tuning (reuse TEM defaults unless overridden)
  const atom_size_mul = opts.atom_size_mul != null ? opts.atom_size_mul : 1.2;
  const atom_size_exp = opts.atom_size_exp != null ? opts.atom_size_exp : 1.25;

  // canvas ctx (optional)
  const canvasCtx = opts.canvasCtx || null;

  const H = img_size[0] | 0;
  const W = img_size[1] | 0;

  const img = newFloatImage(H, W, background_gray);

  if (!atoms || !atoms.length) {
    if (canvasCtx) {
      const imageData = canvasCtx.createImageData(W, H);
      for (let p = 0; p < imageData.data.length; p += 4) {
        imageData.data[p] = background_gray;
        imageData.data[p + 1] = background_gray;
        imageData.data[p + 2] = background_gray;
        imageData.data[p + 3] = 255;
      }
      canvasCtx.putImageData(imageData, 0, 0);
    }
    return new Uint8ClampedArray(H * W);
  }

  const [coords, scale, z_view] = compute_scaled_coordinates(
    atoms,
    img_size,
    angstroms_per_pixel,
    camera,
  );

  const radialCtl = _resolveAFMRadialControls(opts);

  let mergeBonds = [];
  if (bonds == null) {
    const MAX_GUESS_N = 6000;
    if (atoms.length <= MAX_GUESS_N) {
      mergeBonds = guess_bonds_by_distance(atoms);
      try {
        mergeBonds.__tem_guessed = true;
      } catch (_) {}
    }
  } else {
    mergeBonds = bonds;
  }

  // --- 1) Height map (surface topography): max of Gaussian bumps ---
  const height = new Float32Array(H * W);
  const darkField = new Float32Array(H * W);
  let hmax = 0.0;
  let hasDarkOverride = false;

  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    const zv = z_view ? (z_view[i] ?? 0) : (a.z ?? 0);
    if (hide_front && zv > focal_z + 1e-9) continue;

    const x0 = coords[i][0];
    const y0 = coords[i][1];

    const Z = a.Z | 0;
    const rA = get_covalent_radius(Z);

    // Per-element overrides
    const sym = _EO_atomSymbol(a);
    const [sizeMul, darkMul] = _EO_getMul(ov, sym);
    if (Math.abs(darkMul - 1.0) > 1e-9) hasDarkOverride = true;

    // sigma: same geometric logic as TEM atoms (size ~ covalent radius)
    let sigma = Math.max(
      1e-6,
      Math.pow(rA, atom_size_exp) * scale * atom_size_mul,
    );
    sigma *= sizeMul;
    const typical_bond_A = 1.4;
    const cap_rel_bond = 0.55;
    const cap_abs_px = 24.0;
    const cap_px0 = Math.min(
      cap_abs_px,
      Math.max(1.5, cap_rel_bond * typical_bond_A * scale),
    );
    const cap_px = cap_px0 * (Number.isFinite(sizeMul) ? sizeMul : 1.0);
    if (sigma > cap_px) sigma = cap_px;

    const R = Math.ceil(4 * sigma);

    if (x0 + R < 0 || x0 - R >= W || y0 + R < 0 || y0 - R >= H) continue;

    // NOTE:
    //  - height/topography should stay geometric; otherwise dark overrides get swallowed
    //    by global normalization (preview shows ~no change, full image redistributes maxima).
    //  - keep geometry in `height`, and apply darkness later via local contrast field.
    const baseAmp = Math.max(0.05, Math.min(3.0, rA));

    const x0i = Math.floor(x0);
    const y0i = Math.floor(y0);
    const fracX = x0 - x0i;
    const fracY = y0 - y0i;

    const wx = getGaussian1D_cached(sigma, fracX, R);
    const wy = getGaussian1D_cached(sigma, fracY, R);

    const yStart = Math.max(0, y0i - R);
    const yEnd = Math.min(H - 1, y0i + R);
    const xStart = Math.max(0, x0i - R);
    const xEnd = Math.min(W - 1, x0i + R);

    for (let y = yStart; y <= yEnd; y++) {
      const wyv = wy[y - y0i + R];
      const row = y * W;
      for (let x = xStart; x <= xEnd; x++) {
        const g = wx[x - x0i + R] * wyv;
        const bump = baseAmp * g;
        if (bump > 1e-6) {
          const p = row + x;
          if (bump > height[p]) {
            height[p] = bump;
            if (bump > hmax) hmax = bump;
          }
          if (Math.abs(darkMul - 1.0) > 1e-9) {
            darkField[p] += (darkMul - 1.0) * g;
          }
        }
      }
    }
  }

  // normalize height to [0..1]
  if (hmax > 1e-12) {
    const inv = 1.0 / hmax;
    for (let p = 0; p < height.length; p++) height[p] *= inv;
  }

  const mergeData = _buildAFMHeightMergeLayer(W, H, coords, atoms, mergeBonds, z_view, {
    hide_front,
    focal_z,
    atom_size_mul,
    atom_size_exp,
    scale,
  });
  const mergeState = _applyAFMLocalHeightMerge(height, W, H, mergeData);

  const radialField = _buildMoleculeRadialField(
    height,
    W,
    H,
    coords,
    atoms,
    angstroms_per_pixel,
    { profile: radialCtl.profile },
  );

  const globalShapeField = _buildAFMGlobalShapeField(
    height,
    W,
    H,
    mergeData,
    mergeState,
    angstroms_per_pixel,
  );

  const shapeBondTopology =
    mergeBonds && mergeBonds.length
      ? _buildAFMBondTopologyField(W, H, coords, atoms, mergeBonds, z_view, {
          hide_front,
          focal_z,
          atom_size_mul,
          atom_size_exp,
          scale,
          angstroms_per_pixel,
          height,
          mergePresence: mergeState.mergePresence,
          molMask: radialField.molMask,
          shapeSupport: globalShapeField.support,
          shapeContour: globalShapeField.contour,
          shapeOuter: globalShapeField.outer,
          shapeCore: globalShapeField.core,
        })
      : null;

  _applyAFMContourShaping(height, W, H, mergeData, mergeState, shapeBondTopology, radialField, globalShapeField);

  // --- 2) Edge enhancement (simple Laplacian on height) ---
  const edge = new Float32Array(H * W);
  let emax = 0.0;
  if (H >= 3 && W >= 3) {
    for (let y = 1; y < H - 1; y++) {
      const row = y * W;
      for (let x = 1; x < W - 1; x++) {
        const p = row + x;
        const c = height[p];
        const lap =
          4.0 * c -
          height[p - 1] -
          height[p + 1] -
          height[p - W] -
          height[p + W];

        // AFM edge: підсвічуємо тільки «вершинні» області (локальні максимуми) з позитивним Laplacian,
        // щоб уникнути штучних «перемичок» між сусідніми атомами.
        let v = 0.0;
        if (lap > 0.0) {
          const cL = height[p - 1];
          const cR = height[p + 1];
          const cU = height[p - W];
          const cD = height[p + W];
          const isPeak = (c >= cL && c >= cR && c >= cU && c >= cD);
          if (isPeak) {
            const mergedZone = mergeState.mergePresence[p];
            const peakSoft = mergeState.peakSoft[p];
            const edgeKeep = 1.0 - Math.min(0.82, 0.58 * mergedZone + 0.24 * peakSoft);
            v = lap * edgeKeep;
          }
        }
        edge[p] = v;
        if (v > emax) emax = v;
      }
    }
  }
  if (emax > 1e-12) {
    const inv = 1.0 / emax;
    for (let p = 0; p < edge.length; p++) edge[p] *= inv;
  }

  // --- 3) Compose AFM grayscale: bg - K1*height - K2*edge ---
  // Keep coefficients modest; user can tune with Contrast/Background/Invert.
  const K1 = 62.0;
  const K2 = 44.0;
  for (let p = 0; p < img.a.length; p++) {
    img.a[p] = background_gray - (K1 * height[p] + K2 * edge[p]);
  }

  // Local darkness modulation:
  // apply AFTER topography normalization so darkness really changes tone,
  // but BEFORE bonds/postprocess so it does not alter depth slicing semantics.
  if (hasDarkOverride) {
    for (let p = 0; p < img.a.length; p++) {
      let localDark = 1.0 + darkField[p];
      if (localDark < 0.05) localDark = 0.05;
      const d = background_gray - img.a[p];
      img.a[p] = background_gray - d * localDark;
    }
  }

  if (radialCtl.gain > 1e-6) {
    const profile = radialCtl.profile;
    const radialFeatureGain = 0.10 * radialCtl.gain;
    const radialAtomDarkGain = 0.18 * profile.atom_dark_gain * radialCtl.gain;
    const contourProtectArr = shapeBondTopology?.contourField || null;
    const endpointProtectArr = shapeBondTopology?.endpointField || null;
    const openProtectArr = shapeBondTopology?.openField || null;

    for (let p = 0; p < img.a.length; p++) {
      const atomK = radialField.radialAtom[p];
      if (!(atomK > 1e-8)) continue;
      const contourProtect = contourProtectArr ? contourProtectArr[p] : 0.0;
      const endpointProtect = endpointProtectArr ? endpointProtectArr[p] : 0.0;
      const openProtect = openProtectArr ? openProtectArr[p] : 0.0;
      const interiorK = atomK * _clamp01(1.0 - 0.72 * contourProtect - 0.68 * endpointProtect - 0.42 * openProtect);
      if (!(interiorK > 1e-8)) continue;
      const d = background_gray - img.a[p];
      img.a[p] =
        background_gray -
        d * (1.0 + radialAtomDarkGain * interiorK) -
        3.0 * radialFeatureGain * edge[p] * interiorK;
    }

    const radialBlurSigmaPx = profile.blur_sigma_px;
    const radialBlurMix = profile.blur_mix * radialCtl.gain;
    if (radialBlurSigmaPx > 1e-6 && radialBlurMix > 1e-6) {
      const radialBlurred = new Float32Array(img.a);
      const radialScratch = new Float32Array(W * H);
      gaussianBlurFloat(
        { a: radialBlurred, w: W, h: H, bg: 0 },
        radialBlurSigmaPx,
        radialScratch,
      );
      for (let p = 0; p < img.a.length; p++) {
        const mix = radialBlurMix * radialField.radialBlur[p];
        if (!(mix > 1e-8)) continue;
        img.a[p] = img.a[p] * (1.0 - mix) + radialBlurred[p] * mix;
      }
    }
  }

  // --- 4) AFM "bond contrast" overlay (single line per bond; bt ignored) ---
  // Semantics: if bonds are unknown (null), guessing is allowed (with the same safety limits as TEM).
  let use_bonds = [];
  if (draw_bonds_flag) {
    if (bonds == null) {
      const MAX_GUESS_N = 6000;
      const typicalBondA = 1.4;
      const bondPxEst = typicalBondA / Math.max(angstroms_per_pixel, 1e-9);
      const minVisiblePx = Math.max(3.0, bond_wave_width_px * 0.35);
      if (atoms.length <= MAX_GUESS_N && bondPxEst >= minVisiblePx) {
        use_bonds = mergeBonds && mergeBonds.length ? mergeBonds : guess_bonds_by_distance(atoms);
        try {
          use_bonds.__tem_guessed = true;
        } catch (_) {}
      } else {
        use_bonds = [];
      }
    } else {
      use_bonds = bonds;
    }
  }

  // LOD for guessed bonds only (same intent as TEM) — smooth fade, no hard cutoff
  const bondsAreGuessed = bonds == null || !!(bonds && bonds.__tem_guessed);
  let drawBondsEffective = draw_bonds_flag && use_bonds && use_bonds.length > 0;

  // bondK controls how strong bonds are (1..0). Default full strength.
  let bondK = 1.0;

  if (drawBondsEffective && bondsAreGuessed) {
    const maxSamples = 256;
    const step = Math.max(1, Math.floor(use_bonds.length / maxSamples));
    let sum = 0,
      cnt = 0;

    for (let k = 0; k < use_bonds.length; k += step) {
      const b = use_bonds[k];
      const i = b[0],
        j = b[1];
      const p1 = coords[i],
        p2 = coords[j];
      if (!p1 || !p2) continue;

      const dx = p2[0] - p1[0];
      const dy = p2[1] - p1[1];
      const L = Math.hypot(dx, dy);
      if (Number.isFinite(L)) {
        sum += L;
        cnt++;
      }

      if (cnt >= maxSamples) break;
    }

    const avgLenPx = cnt ? sum / cnt : 0;
    const minVisiblePx = Math.max(3.0, bond_wave_width_px * 0.35);

    // Fade window: above minVisiblePx -> full; below ~0.55*minVisiblePx -> off
    const fadeStart = minVisiblePx;
    const fadeEnd = minVisiblePx * 0.55;

    if (avgLenPx <= fadeEnd) {
      bondK = 0.0;
    } else if (avgLenPx >= fadeStart) {
      bondK = 1.0;
    } else {
      bondK = (avgLenPx - fadeEnd) / (fadeStart - fadeEnd); // linear fade
    }

    // If essentially invisible, skip drawing for performance
    if (bondK <= 1e-3) {
      drawBondsEffective = false;
      use_bonds = [];
      bondK = 0.0;
    }
  }

  if (!drawBondsEffective) {
    use_bonds = [];
    bondK = 0.0;
  }

  if (drawBondsEffective && use_bonds && use_bonds.length) {
    const bondLayer = { a: new Float32Array(H * W), h: H, w: W, bg: 0 };

    draw_bonds_single(bondLayer, coords, use_bonds, atoms, {
      wave_width_px: bond_wave_width_px,
      wave_amp: bond_wave_amplitude,
      focal_z,
      hide_front,
      z_view,
    });

    const profile = radialCtl.profile;
    const bondTopology =
      shapeBondTopology && use_bonds === mergeBonds
        ? shapeBondTopology
        : _buildAFMBondTopologyField(W, H, coords, atoms, use_bonds, z_view, {
            hide_front,
            focal_z,
            atom_size_mul,
            atom_size_exp,
            scale,
            angstroms_per_pixel,
            height,
            mergePresence: mergeState.mergePresence,
            molMask: radialField.molMask,
            shapeSupport: globalShapeField.support,
            shapeContour: globalShapeField.contour,
            shapeOuter: globalShapeField.outer,
            shapeCore: globalShapeField.core,
          });

    const bondComp = _buildAFMBondBlurCompLayers(
      bondLayer,
      W,
      H,
      blur_sigma,
      bond_wave_width_px,
      {
        height,
        mergePresence: mergeState.mergePresence,
        molMask: radialField.molMask,
        bondFeature: bondTopology.feature,
        bondMask: bondTopology.mask,
        openField: bondTopology.openField,
        contourField: bondTopology.contourField,
        endpointField: bondTopology.endpointField,
      },
    );

    for (let p = 0; p < bondLayer.a.length; p++) {
      const bondBase = _composeAFMBondValue(p, bondComp, {
        height,
        mergePresence: mergeState.mergePresence,
        molMask: radialField.molMask,
        bondFeature: bondTopology.feature,
        bondMask: bondTopology.mask,
      });
      if (!(bondBase > 1e-8)) continue;

      const topoK = bondTopology.feature[p];
      const tipK = bondTopology.tipField[p];
      const bendK = bondTopology.bendField[p];
      const endpointFieldK = bondTopology.endpointField ? bondTopology.endpointField[p] : 0.0;
      const cornerFieldK = bondTopology.cornerField ? bondTopology.cornerField[p] : 0.0;
      const branchFieldK = bondTopology.branchField ? bondTopology.branchField[p] : 0.0;
      const terminalFieldK = bondTopology.terminalField ? bondTopology.terminalField[p] : 0.0;
      const terminalSpillFieldK = bondTopology.terminalSpillField ? bondTopology.terminalSpillField[p] : 0.0;
      const openFieldK = bondTopology.openField ? bondTopology.openField[p] : 0.0;
      const contourFieldK = bondTopology.contourField ? bondTopology.contourField[p] : 0.0;

      const terminalK = _smoothstep(0.12, 0.52, terminalFieldK);
      const terminalSpillK = _smoothstep(0.08, 0.44, terminalSpillFieldK);
      const terminalEnvK = _clamp01(Math.max(terminalK, 0.70 * terminalSpillK));
      const selectiveK = _smoothstep(
        0.14,
        0.50,
        Math.max(0.42 * topoK + 0.58 * terminalEnvK, 0.48 * openFieldK + 0.42 * contourFieldK),
      ) * (0.66 + 0.34 * terminalEnvK);
      const endpointK = _smoothstep(0.12, 0.48, Math.max(endpointFieldK, 0.46 * tipK, 0.46 * bendK)) * (0.72 + 0.28 * terminalEnvK);
      const cornerK = _smoothstep(0.12, 0.50, Math.max(cornerFieldK, 0.54 * bendK)) * (0.70 + 0.30 * terminalEnvK);
      const openK = _smoothstep(0.14, 0.52, openFieldK) * (0.68 + 0.32 * terminalEnvK);
      const branchK = _smoothstep(0.18, 0.68, branchFieldK) * (0.74 + 0.26 * terminalEnvK);
      const denseDamp =
        1.0 -
        0.18 *
          _smoothstep(0.62, 0.96, Math.max(height[p], radialField.molMask[p])) *
          (1.0 - 0.62 * contourFieldK);

      const contourK = _smoothstep(0.14, 0.52, contourFieldK) * (0.66 + 0.34 * terminalEnvK);
      const baseGain = 0.28 + 0.36 * Math.max(0.0, profile.base_bond_gain);
      const selectiveGain = 0.70 * profile.bond_gain * selectiveK;
      const endpointGain = 0.56 * profile.bond_gain * endpointK;
      const cornerGain = 0.32 * profile.bond_gain * cornerK;
      const openGain = 0.36 * profile.bond_gain * openK * (1.0 - 0.10 * branchK);
      const contourGain = 0.24 * profile.bond_gain * contourK;
      const spillGain = 0.22 * profile.bond_gain * terminalSpillK * (1.0 - 0.58 * terminalK);
      const manualGain = 0.18 * Math.max(0.0, radialCtl.bondGain);
      const blurComp =
        1.0 +
        0.10 * bondComp.blurK +
        0.12 * bondComp.blurK * Math.max(selectiveK, openK, contourK, 0.7 * endpointK, 0.6 * terminalSpillK);
      const debugDetectBoost = 1.0;
      const preBoost =
        (baseGain + selectiveGain + endpointGain + cornerGain + openGain + contourGain + spillGain + manualGain) *
        blurComp *
        denseDamp *
        (0.82 + 0.18 * terminalEnvK) *
        debugDetectBoost;

      img.a[p] += bondBase * bondK * preBoost;
    }
  }

  // --- 5) Postprocess: blur/noise/contrast/invert/clipping (same knobs as TEM; no DoF here) ---
  const scratch = new Float32Array(H * W);
  if (blur_sigma > 0) gaussianBlurFloat(img, blur_sigma, scratch);

  if (noise_stddev > 0) {
    const σ = noise_stddev;
    for (let i = 0; i < img.a.length; i++) {
      const u1 = Math.random(),
        u2 = Math.random();
      const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
      img.a[i] += z * σ;
    }
  }

  if (contrast !== 1.0) {
    for (let i = 0; i < img.a.length; i++) {
      img.a[i] = (img.a[i] - background_gray) * contrast + background_gray;
    }
  }

  if (invert) {
    for (let i = 0; i < img.a.length; i++) {
      const d = background_gray - img.a[i];
      img.a[i] = background_gray + d;
    }
  }

  if (low_clip != null || high_clip != null) {
    const lo = (low_clip == null ? 0 : low_clip) | 0;
    const hi = (high_clip == null ? 255 : high_clip) | 0;
    const a = Math.min(lo, hi),
      b = Math.max(lo, hi);
    for (let i = 0; i < img.a.length; i++) {
      const v = img.a[i];
      img.a[i] = v < a ? a : v > b ? b : v;
    }
  }

  // to 8-bit
  const out = new Uint8ClampedArray(H * W);
  for (let i = 0; i < out.length; i++) out[i] = clamp255(img.a[i]);

  // draw to canvas
  if (canvasCtx) {
    const imageData = canvasCtx.createImageData(W, H);
    for (let i = 0, p = 0; i < out.length; i++, p += 4) {
      const v = out[i];
      imageData.data[p] = v;
      imageData.data[p + 1] = v;
      imageData.data[p + 2] = v;
      imageData.data[p + 3] = 255;
    }
    canvasCtx.putImageData(imageData, 0, 0);

    if (show_scale_bar) {
      draw_scale_bar(canvasCtx, { h: H, w: W }, angstroms_per_pixel, {
        corner: scale_bar_corner,
        margin: scale_bar_margin_px,
        invert,
      });
    }
  }

  return out;
}
