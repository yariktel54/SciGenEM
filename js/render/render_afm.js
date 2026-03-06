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
  let s = a.sym || a.el || a.symbol || a.element;
  if (typeof s === "string") {
    s = s.trim();
    if (s) return s;
  }
  const Z = a.Z;
  if (Number.isFinite(Z)) {
    const zi = Z | 0;
    if (zi > 0 && zi < _EO_Z2SYM.length) return _EO_Z2SYM[zi] || null;
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

  // --- 1) Height map (surface topography): max of Gaussian bumps ---
  const height = new Float32Array(H * W);
  let hmax = 0.0;

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

    // sigma: same geometric logic as TEM atoms (size ~ covalent radius)
    let sigma = Math.max(
      1e-6,
      Math.pow(rA, atom_size_exp) * scale * atom_size_mul,
    );
    sigma *= sizeMul;
    const typical_bond_A = 1.4;
    const cap_rel_bond = 0.55;
    const cap_abs_px = 24.0;
    const cap_px = Math.min(
      cap_abs_px,
      Math.max(1.5, cap_rel_bond * typical_bond_A * scale),
    );
    if (sigma > cap_px) sigma = cap_px;

    const R = Math.ceil(4 * sigma);

    if (x0 + R < 0 || x0 - R >= W || y0 + R < 0 || y0 - R >= H) continue;

    const amp = darkMul * Math.max(0.05, Math.min(3.0, rA)); // size-driven; later normalized by hmax

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
        const bump = amp * (wx[x - x0i + R] * wyv);
        if (bump <= 1e-6) continue;
        const p = row + x;
        if (bump > height[p]) {
          height[p] = bump;
          if (bump > hmax) hmax = bump;
        }
      }
    }
  }

  // normalize height to [0..1]
  if (hmax > 1e-12) {
    const inv = 1.0 / hmax;
    for (let p = 0; p < height.length; p++) height[p] *= inv;
  }

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
            v = lap;
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
        use_bonds = guess_bonds_by_distance(atoms);
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

    // add to AFM base as a bright "bond contrast" layer; global invert stays global
    for (let p = 0; p < bondLayer.a.length; p++) {
      img.a[p] += bondLayer.a[p] * bondK;
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
