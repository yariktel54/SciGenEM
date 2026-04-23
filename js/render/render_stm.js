// js/render/render_stm.js
// STM-like surface render rebuilt on top of AFM-style envelope logic.
//
// Design intent:
//  - STM must stay surface-sensitive (not volume / x-ray-like)
//  - apparent contrast should follow the top accessible electronic envelope,
//    not a sum through the full thickness
//  - atoms should appear as localized bright protrusions / lobes on the surface,
//    not as literal TEM projected columns or dark pits
//  - bonds are NOT drawn in STM mode
//  - DoF is ignored in STM mode; hide_front/focal_z may still act as an educational slice

import {
    compute_scaled_coordinates,
    newFloatImage,
    gaussianBlurFloat,
    clamp255,
    getGaussian1D_cached,
    draw_scale_bar,
    draw_stm_orange_frame_to_canvas
} from './render_common.js';

// ---------------------------
// Helpers
// ---------------------------

let _randn_spare = null;
function randn() {
    if (_randn_spare != null) {
        const v = _randn_spare;
        _randn_spare = null;
        return v;
    }
    let u1 = 0, u2 = 0;
    while (u1 === 0) u1 = Math.random();
    u2 = Math.random();
    const r = Math.sqrt(-2.0 * Math.log(u1));
    const t = 2.0 * Math.PI * u2;
    _randn_spare = r * Math.sin(t);
    return r * Math.cos(t);
}

function meanStd(arr) {
    let m = 0;
    for (let i = 0; i < arr.length; i++) m += arr[i];
    m /= Math.max(1, arr.length);
    let v = 0;
    for (let i = 0; i < arr.length; i++) {
        const d = arr[i] - m;
        v += d * d;
    }
    v /= Math.max(1, arr.length);
    return [m, Math.sqrt(Math.max(0, v))];
}

function normalizeStdInPlace(arr, eps = 1e-9) {
    const [m, s] = meanStd(arr);
    const inv = 1.0 / Math.max(eps, s);
    for (let i = 0; i < arr.length; i++) arr[i] = (arr[i] - m) * inv;
    return arr;
}

function gaussianSmooth1D_inplace(arr, sigma) {
    if (!(sigma > 0)) return arr;
    const n = arr.length;
    const ksize = Math.max(3, (Math.floor(sigma * 3) * 2 + 1));
    const half = ksize >> 1;

    const kernel = new Float32Array(ksize);
    const inv2sig2 = 1.0 / (2.0 * sigma * sigma);
    let sum = 0;
    for (let k = -half; k <= half; k++) {
        const v = Math.exp(-(k * k) * inv2sig2);
        kernel[k + half] = v;
        sum += v;
    }
    for (let i = 0; i < ksize; i++) kernel[i] /= sum;

    const tmp = new Float32Array(n);
    for (let i = 0; i < n; i++) {
        let s = 0;
        for (let k = -half; k <= half; k++) {
            const j = (i + k < 0) ? 0 : (i + k >= n) ? (n - 1) : (i + k);
            s += arr[j] * kernel[k + half];
        }
        tmp[i] = s;
    }
    arr.set(tmp);
    return arr;
}

function randomWalk1D(n, stepSigma = 1.0) {
    const out = new Float32Array(n);
    let v = 0;
    for (let i = 0; i < n; i++) {
        v += randn() * stepSigma;
        out[i] = v;
    }
    return out;
}

function shiftRowLinear(src, srcOff, dst, w, dx) {
    if (!dx) {
        for (let x = 0; x < w; x++) dst[x] = src[srcOff + x];
        return;
    }
    for (let x = 0; x < w; x++) {
        const sx = x - dx;
        if (sx <= 0) { dst[x] = src[srcOff + 0]; continue; }
        if (sx >= (w - 1)) { dst[x] = src[srcOff + (w - 1)]; continue; }
        const x0 = sx | 0;
        const t = sx - x0;
        const v0 = src[srcOff + x0];
        const v1 = src[srcOff + (x0 + 1)];
        dst[x] = v0 + (v1 - v0) * t;
    }
}

function buildPeakOnlyEdge(height, H, W) {
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

                let v = 0.0;
                if (lap > 0.0) {
                    const cL = height[p - 1];
                    const cR = height[p + 1];
                    const cU = height[p - W];
                    const cD = height[p + W];
                    const isPeak = (c >= cL && c >= cR && c >= cU && c >= cD);
                    if (isPeak) v = lap;
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
    return edge;
}

function normalizeField(field) {
    let fmax = 0.0;
    for (let i = 0; i < field.length; i++) if (field[i] > fmax) fmax = field[i];
    if (fmax > 1e-12) {
        const inv = 1.0 / fmax;
        for (let i = 0; i < field.length; i++) field[i] *= inv;
    }
    return field;
}

function applyScanlineArtifacts(sigA, H, W, opts) {
    const scan_off_amp = (opts.stm_scan_off_amp != null) ? Number(opts.stm_scan_off_amp) : 0.018;
    const scan_gain_amp = (opts.stm_scan_gain_amp != null) ? Number(opts.stm_scan_gain_amp) : 0.055;
    const scan_band_amp = (opts.stm_scan_band_amp != null) ? Number(opts.stm_scan_band_amp) : 0.015;
    const scan_drift_amp_px = (opts.stm_scan_drift_amp_px != null) ? Number(opts.stm_scan_drift_amp_px) : 0.95;
    const scan_smooth_sigma = (opts.stm_scan_smooth_sigma != null) ? Number(opts.stm_scan_smooth_sigma) : 10.0;

    let off = randomWalk1D(H, 1.0);
    let gain = randomWalk1D(H, 1.0);
    let dx = randomWalk1D(H, 1.0);

    gaussianSmooth1D_inplace(off, scan_smooth_sigma);
    gaussianSmooth1D_inplace(gain, scan_smooth_sigma);
    gaussianSmooth1D_inplace(dx, scan_smooth_sigma);

    normalizeStdInPlace(off);
    normalizeStdInPlace(gain);
    normalizeStdInPlace(dx);

    for (let y = 0; y < H; y++) {
        off[y] *= scan_off_amp;
        gain[y] *= scan_gain_amp;
        dx[y] *= scan_drift_amp_px;
    }

    const per1 = Math.max(24, Math.round(H * 0.42));
    const per2 = Math.max(36, Math.round(H * 0.73));
    const phi = Math.random() * Math.PI * 2;

    const rowBuf = new Float32Array(W);
    for (let y = 0; y < H; y++) {
        const band = scan_band_amp * (
            Math.sin((2 * Math.PI * y) / per1 + phi) +
            0.45 * Math.sin((2 * Math.PI * y) / per2 + 1.7 * phi)
        );

        const rowOff = y * W;
        shiftRowLinear(sigA, rowOff, rowBuf, W, dx[y]);

        const g = 1.0 + gain[y];
        const o = off[y] + band;
        for (let x = 0; x < W; x++) {
            let v = rowBuf[x] * g + o;
            if (v < 0) v = 0;
            sigA[rowOff + x] = v;
        }
    }

    normalizeField(sigA);
}

function getScanlinesFlag(opts) {
    if (!opts) return false;
    const direct = (opts.cb_scanlines ?? opts.cbScanlines ?? opts.scanlines_flag ?? opts.scanlines ?? opts.enable_scanlines);
    if (typeof direct === 'boolean') return direct;
    const ui = opts.ui || opts.controls || opts.flags || null;
    if (ui) {
        const v = (ui.cb_scanlines ?? ui.cbScanlines ?? ui.scanlines_flag ?? ui.scanlines ?? ui.enable_scanlines);
        if (typeof v === 'boolean') return v;
    }
    return false;
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

function estimateStmGeom(atom, scale, ov, opts) {
    const Z = Math.max(1, atom?.Z | 0);
    const sym = _EO_atomSymbol(atom);
    const [sizeMul, darkMul] = _EO_getMul(ov, sym);

    // STM apparent lateral extent is primarily tip / vacuum / electronic-envelope limited,
    // not a universal monotonic function of atomic number.
    // We therefore use a global tip-limited envelope width and keep element overrides manual only.
    const tipRadiusA = Math.max(0.20, Number(opts.stm_tip_radius_A ?? 3.0));
    const gapA = Math.max(0.20, Number(opts.stm_gap_A ?? 3.5));
    const decayLenA = Math.max(0.20, Number(opts.stm_decay_length_A ?? 1.0));
    const sigmaFromResolution = Math.max(0.05, Number(opts.stm_sigma_from_resolution ?? 0.22));
    const sizeOverrideStrength = Math.max(0.0, Number(opts.stm_size_override_strength ?? 1.0));
    const typicalBondA = Math.max(0.25, Number(opts.typical_bond_A ?? 1.40));
    const capRelBond = Math.max(0.20, Number(opts.stm_sigma_cap_rel_bond ?? 0.85));
    const capAbsPx = Math.max(1.0, Number(opts.stm_sigma_cap_abs_px ?? 28.0));

    // Tersoff-Hamann lateral resolution scale ~ sqrt[(2 A) * (R + d)].
    // Here we compress that scale into a Gaussian sigma for a compact educational blob model.
    const tipLimitedResolutionA = Math.sqrt(Math.max(1e-6, 2.0 * decayLenA * (tipRadiusA + gapA)));
    let sigmaA = Number(opts.stm_sigma_A);
    if (!(sigmaA > 0)) sigmaA = tipLimitedResolutionA * sigmaFromResolution;

    let sigma = Math.max(1e-6, sigmaA * scale);

    // Keep explicit per-element size tweaks available only as manual artistic overrides.
    const manualSizeMul = 1.0 + (sizeMul - 1.0) * sizeOverrideStrength;
    sigma *= manualSizeMul;
    if (sigma < 0.45) sigma = 0.45;

    const capPx = Math.min(capAbsPx, Math.max(1.20, capRelBond * typicalBondA * scale));
    if (sigma > capPx) sigma = capPx;

    return { sigma, sizeMul, darkMul, tipLimitedResolutionA, sigmaA, Z };
}

function buildSurfaceSubset(atoms, coords, z_view, scale, ov, opts, H, W) {
    const cellPx = Math.max(2.0, Number(opts.stm_surface_cell_px ?? 4.0));
    const gridW = Math.max(1, Math.ceil(W / cellPx));
    const gridH = Math.max(1, Math.ceil(H / cellPx));
    const topZ = new Float32Array(gridW * gridH);
    const NEG = -1e30;
    for (let i = 0; i < topZ.length; i++) topZ[i] = NEG;

    const footprintMul = Math.max(0.75, Number(opts.stm_surface_footprint_mul ?? 1.20));
    const envelopeA = Math.max(0.0, Number(opts.stm_surface_envelope_A ?? 0.70));
    const padPx = 0.5 * cellPx;
    const candidates = [];

    function cellCoverLoop(x0, y0, radiusPx, fn) {
        const gx0 = Math.max(0, Math.floor((x0 - radiusPx) / cellPx));
        const gx1 = Math.min(gridW - 1, Math.floor((x0 + radiusPx) / cellPx));
        const gy0 = Math.max(0, Math.floor((y0 - radiusPx) / cellPx));
        const gy1 = Math.min(gridH - 1, Math.floor((y0 + radiusPx) / cellPx));
        const rr = radiusPx + padPx;
        const rr2 = rr * rr;

        for (let gy = gy0; gy <= gy1; gy++) {
            const cy = (gy + 0.5) * cellPx;
            for (let gx = gx0; gx <= gx1; gx++) {
                const cx = (gx + 0.5) * cellPx;
                const dx = cx - x0;
                const dy = cy - y0;
                if ((dx * dx + dy * dy) > rr2) continue;
                fn(gx, gy);
            }
        }
    }

    for (let i = 0; i < atoms.length; i++) {
        const p = coords[i];
        if (!p) continue;
        const x0 = p[0];
        const y0 = p[1];
        const zAtom = z_view ? (z_view[i] ?? 0) : (atoms[i].z ?? 0);
        if (opts.hide_front && zAtom > opts.focal_z + 1e-9) continue;

        const geom = estimateStmGeom(atoms[i], scale, ov, opts);
        const surfaceRadiusPx = Math.max(cellPx * 0.8, geom.sigma * footprintMul);
        if ((x0 + surfaceRadiusPx) < 0 || (x0 - surfaceRadiusPx) >= W || (y0 + surfaceRadiusPx) < 0 || (y0 - surfaceRadiusPx) >= H) continue;

        const cand = {
            index: i,
            x0,
            y0,
            z: zAtom,
            sigma: geom.sigma,
            sizeMul: geom.sizeMul,
            darkMul: geom.darkMul,
            Z: geom.Z,
            surfaceRadiusPx,
            localTopZ: zAtom,
            surfaceGapA: Infinity,
            selected: false,
        };
        candidates.push(cand);

        cellCoverLoop(x0, y0, surfaceRadiusPx, function (gx, gy) {
            const gi = gy * gridW + gx;
            if (zAtom > topZ[gi]) topZ[gi] = zAtom;
        });
    }

    if (!candidates.length) return [];

    for (let i = 0; i < candidates.length; i++) {
        const cand = candidates[i];
        let localTop = NEG;
        let minGap = Infinity;
        cellCoverLoop(cand.x0, cand.y0, cand.surfaceRadiusPx, function (gx, gy) {
            const zTop = topZ[gy * gridW + gx];
            if (zTop <= NEG * 0.5) return;
            if (zTop > localTop) localTop = zTop;
            const gap = zTop - cand.z;
            if (gap < minGap) minGap = gap;
            if (gap <= envelopeA + 1e-9) cand.selected = true;
        });

        if (localTop > NEG * 0.5) cand.localTopZ = localTop;
        if (Number.isFinite(minGap)) cand.surfaceGapA = Math.max(0.0, minGap);
    }

    let anySelected = false;
    for (let i = 0; i < candidates.length; i++) {
        if (candidates[i].selected) {
            anySelected = true;
            break;
        }
    }
    if (!anySelected) {
        let best = 0;
        for (let i = 1; i < candidates.length; i++) {
            if (candidates[i].z > candidates[best].z) best = i;
        }
        candidates[best].selected = true;
        candidates[best].surfaceGapA = 0.0;
        candidates[best].localTopZ = candidates[best].z;
    }

    const out = [];
    for (let i = 0; i < candidates.length; i++) {
        if (candidates[i].selected) out.push(candidates[i]);
    }
    return out;
}

function depositMaxGaussian(field, H, W, x0, y0, sigma, amp) {
    if (!(sigma > 1e-9) || !(amp > 0)) return;

    const x0i = Math.floor(x0);
    const y0i = Math.floor(y0);
    const fracX = x0 - x0i;
    const fracY = y0 - y0i;
    const R = Math.ceil(4 * sigma);

    if ((x0 + R) < 0 || (x0 - R) >= W || (y0 + R) < 0 || (y0 - R) >= H) return;

    const wx = getGaussian1D_cached(sigma, fracX, R);
    const wy = getGaussian1D_cached(sigma, fracY, R);

    const yStart = Math.max(0, y0i - R);
    const yEnd = Math.min(H - 1, y0i + R);
    const xStart = Math.max(0, x0i - R);
    const xEnd = Math.min(W - 1, x0i + R);

    for (let y = yStart; y <= yEnd; y++) {
        const wyv = wy[(y - y0i) + R];
        const row = y * W;
        for (let x = xStart; x <= xEnd; x++) {
            const v = amp * (wx[(x - x0i) + R] * wyv);
            const p = row + x;
            if (v > field[p]) field[p] = v;
        }
    }
}

function accumulateGaussian(field, H, W, x0, y0, sigma, amp) {
    if (!(sigma > 1e-9) || !(amp !== 0)) return;

    const x0i = Math.floor(x0);
    const y0i = Math.floor(y0);
    const fracX = x0 - x0i;
    const fracY = y0 - y0i;
    const R = Math.ceil(4 * sigma);

    if ((x0 + R) < 0 || (x0 - R) >= W || (y0 + R) < 0 || (y0 - R) >= H) return;

    const wx = getGaussian1D_cached(sigma, fracX, R);
    const wy = getGaussian1D_cached(sigma, fracY, R);

    const yStart = Math.max(0, y0i - R);
    const yEnd = Math.min(H - 1, y0i + R);
    const xStart = Math.max(0, x0i - R);
    const xEnd = Math.min(W - 1, x0i + R);

    for (let y = yStart; y <= yEnd; y++) {
        const wyv = wy[(y - y0i) + R];
        const row = y * W;
        for (let x = xStart; x <= xEnd; x++) {
            field[row + x] += amp * (wx[(x - x0i) + R] * wyv);
        }
    }
}

export function render_stm_like(atoms, opts = {}) {
    const ov = opts?.element_overrides || null;
    const {
        img_size = [400, 400],
        angstroms_per_pixel = 0.1,
        blur_sigma = 1.0,
        background_gray = 127,
        invert = false,
        noise_stddev = 2.0,
        contrast = 1.0,
        camera = null,
        low_clip = null,
        high_clip = null,
        focal_z = 0.0,
        hide_front = false,
        show_scale_bar = false,
        scale_bar_corner = 'bl',
        scale_bar_margin_px = 12,
        canvasCtx = null
    } = opts;

    const [H, W] = img_size;
    if (!atoms || !atoms.length) {
        const out0 = new Uint8ClampedArray(H * W);
        out0.fill(clamp255(background_gray));
        if (canvasCtx) {
            draw_stm_orange_frame_to_canvas(canvasCtx, out0, W, H);
            if (show_scale_bar) {
                draw_scale_bar(canvasCtx, { h: H, w: W }, angstroms_per_pixel, {
                    corner: scale_bar_corner,
                    margin: scale_bar_margin_px,
                    invert
                });
            }
        }
        return out0;
    }

    const [coords, scale, z_view] = compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel, camera);
    const surfaceAtoms = buildSurfaceSubset(atoms, coords, z_view, scale, ov, {
        ...opts,
        focal_z,
        hide_front,
    }, H, W);

    const height = new Float32Array(H * W);
    const center = new Float32Array(H * W);
    const darkField = new Float32Array(H * W);
    let hasDarkOverride = false;

    let topSurfaceZ = -Infinity;
    for (let i = 0; i < surfaceAtoms.length; i++) {
        if (surfaceAtoms[i].localTopZ > topSurfaceZ) topSurfaceZ = surfaceAtoms[i].localTopZ;
    }
    if (!Number.isFinite(topSurfaceZ)) topSurfaceZ = 0.0;

    const stm_surface_alpha = (opts.stm_surface_alpha != null) ? Number(opts.stm_surface_alpha) : 3.6;
    const stm_global_alpha = (opts.stm_global_alpha != null) ? Number(opts.stm_global_alpha) : 0.28;
    const stm_surface_floor = (opts.stm_surface_floor != null) ? Number(opts.stm_surface_floor) : 0.55;
    const stm_terrace_floor = (opts.stm_terrace_floor != null) ? Number(opts.stm_terrace_floor) : 0.45;
    const stm_center_sigma_mul = (opts.stm_center_sigma_mul != null) ? Number(opts.stm_center_sigma_mul) : 0.52;

    for (let i = 0; i < surfaceAtoms.length; i++) {
        const s = surfaceAtoms[i];
        const zGapLocal = Math.max(0.0, s.surfaceGapA);
        const zGapGlobal = Math.max(0.0, topSurfaceZ - s.z);

        const wSurfaceRaw = Math.exp(-stm_surface_alpha * zGapLocal);
        const wTerraceRaw = Math.exp(-stm_global_alpha * zGapGlobal);
        const wSurface = stm_surface_floor + (1.0 - stm_surface_floor) * wSurfaceRaw;
        const wTerrace = stm_terrace_floor + (1.0 - stm_terrace_floor) * wTerraceRaw;
        const w = wSurface * wTerrace;
        if (w < 1e-5) continue;

        // STM apparent contrast should not follow a hard-coded monotonic atomic-number law.
        // In the simplified surface model we keep visibility mostly geometric/surface-driven;
        // per-element contrast can still be adjusted via element overrides if needed.
        const geomAmp = 0.82 * w;
        const centerAmp = 0.98 * w;

        depositMaxGaussian(height, H, W, s.x0, s.y0, s.sigma, geomAmp);
        depositMaxGaussian(center, H, W, s.x0, s.y0, Math.max(0.35, s.sigma * stm_center_sigma_mul), centerAmp);

        if (Math.abs(s.darkMul - 1.0) > 1e-9) {
            hasDarkOverride = true;
            accumulateGaussian(darkField, H, W, s.x0, s.y0, s.sigma, (s.darkMul - 1.0));
        }
    }

    normalizeField(height);
    normalizeField(center);
    const edge = buildPeakOnlyEdge(height, H, W);

    // STM apparent contrast: broad surface envelope + localized electronic apex.
    // Kept on AFM-like max/envelope semantics, but without AFM bond/topography styling.
    const sig = newFloatImage(H, W, 0);
    const sigA = sig.a;
    const Kheight = (opts.stm_k_height != null) ? Number(opts.stm_k_height) : 0.72;
    const Kcenter = (opts.stm_k_center != null) ? Number(opts.stm_k_center) : 0.68;
    const Kedge = (opts.stm_k_edge != null) ? Number(opts.stm_k_edge) : 0.16;

    for (let p = 0; p < sigA.length; p++) {
        sigA[p] = Kheight * height[p] + Kcenter * center[p] + Kedge * edge[p];
    }
    normalizeField(sigA);

    const stm_mid_gamma = (opts.stm_mid_gamma != null) ? Number(opts.stm_mid_gamma) : 0.72;
    if (stm_mid_gamma > 0 && Math.abs(stm_mid_gamma - 1.0) > 1e-9) {
        for (let p = 0; p < sigA.length; p++) {
            const v = sigA[p];
            sigA[p] = (v > 0) ? Math.pow(v, stm_mid_gamma) : 0;
        }
        normalizeField(sigA);
    }

    if (getScanlinesFlag(opts)) {
        applyScanlineArtifacts(sigA, H, W, opts);
    }

    // Constant-current STM topography is closer to a logarithmic apparent-height map
    // than to a darkened absorption image. Convert the surface signal into a bright
    // protrusion field with gentle logarithmic compression.
    const stm_topo_log_gain = (opts.stm_topo_log_gain != null) ? Number(opts.stm_topo_log_gain) : 5.5;
    const topoField = new Float32Array(sigA.length);
    const topoNorm = Math.log1p(Math.max(1e-6, stm_topo_log_gain));
    for (let p = 0; p < sigA.length; p++) {
        const v = Math.max(0.0, sigA[p]);
        topoField[p] = Math.log1p(stm_topo_log_gain * v) / topoNorm;
    }
    normalizeField(topoField);

    const stm_signal_strength = (opts.stm_signal_strength != null) ? Number(opts.stm_signal_strength) : 228.0;
    const img = newFloatImage(H, W, background_gray);
    for (let p = 0; p < img.a.length; p++) {
        img.a[p] = background_gray + (stm_signal_strength * topoField[p]);
    }

    if (hasDarkOverride) {
        const stm_dark_override_strength = (opts.stm_dark_override_strength != null) ? Number(opts.stm_dark_override_strength) : 0.18;
        for (let p = 0; p < img.a.length; p++) {
            let localMul = 1.0 - stm_dark_override_strength * darkField[p];
            if (localMul < 0.20) localMul = 0.20;
            img.a[p] = background_gray + (img.a[p] - background_gray) * localMul;
        }
    }

    const scratch = new Float32Array(H * W);
    if (blur_sigma > 0) gaussianBlurFloat(img, blur_sigma, scratch);

    if (noise_stddev > 0) {
        const sigmaNoise = noise_stddev;
        for (let i = 0; i < img.a.length; i++) img.a[i] += randn() * sigmaNoise;
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
        const a = Math.min(lo, hi), b = Math.max(lo, hi);
        for (let i = 0; i < img.a.length; i++) {
            const v = img.a[i];
            img.a[i] = v < a ? a : v > b ? b : v;
        }
    }

    const out = new Uint8ClampedArray(H * W);
    for (let i = 0; i < out.length; i++) out[i] = clamp255(img.a[i]);

    if (canvasCtx) {
        draw_stm_orange_frame_to_canvas(canvasCtx, out, W, H);

        if (show_scale_bar) {
            draw_scale_bar(canvasCtx, { h: H, w: W }, angstroms_per_pixel, {
                corner: scale_bar_corner,
                margin: scale_bar_margin_px,
                invert
            });
        }
    }

    return out;
}

export function render_image_stm(atoms, opts = {}) {
    return render_stm_like(atoms, opts);
}
