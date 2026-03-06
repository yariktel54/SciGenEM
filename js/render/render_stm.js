// js/render/render_stm.js
// Wave 3: realistic STM-like visualization with true scanlines/drift
//
// Rules:
//  - bonds are NOT drawn in STM mode
//  - DoF is ignored (surface mode)
//  - hide_front + focal_z act as a "slice": if hide_front==true, atoms with z > focal_z are excluded
//
// Implementation notes:
//  - Base signal is a fast "current map": sum of 2D Gaussians with strong exponential falloff by z
//  - Optional scanline block adds: per-row DC offset, per-row gain, low-frequency banding, and gentle x-drift per row
//  - Postprocess matches TEM/AFM pipeline: blur/noise/contrast/invert/clipping + grayscale->RGBA

import { get_covalent_radius } from '../bonds/bonds.js';
import {
    compute_scaled_coordinates,
    newFloatImage,
    gaussianBlurFloat,
    clamp255,
    getGaussian1D_cached,
    draw_scale_bar
} from './render_common.js';

// ---------------------------
// Helpers
// ---------------------------

// Fast randn (Box–Muller with cache)
let _randn_spare = null;
function randn() {
    if (_randn_spare != null) {
        const v = _randn_spare;
        _randn_spare = null;
        return v;
    }
    let u1 = 0, u2 = 0;
    // avoid log(0)
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

// Shift one row by dx (subpixel) with linear interpolation.
// dx > 0 shifts features to the right.
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

// Add a positive 2D Gaussian bump into sig (Float32Array) using cached 1D weights.
// (ROI-based, so O(N * R^2), not O(N*W*H).)
function addGaussianROI(sigA, H, W, x0, y0, sigma, amp) {
    if (!(sigma > 1e-9) || amp === 0) return;

    const x0i = Math.floor(x0);
    const y0i = Math.floor(y0);
    const fracX = x0 - x0i;
    const fracY = y0 - y0i;

    const R = Math.ceil(5 * sigma);

    // fully outside
    if ((x0 + R) < 0 || (x0 - R) >= W || (y0 + R) < 0 || (y0 - R) >= H) return;

    const wx = getGaussian1D_cached(sigma, fracX, R);
    const wy = getGaussian1D_cached(sigma, fracY, R);

    const yStart = Math.max(0, (y0i - R));
    const yEnd   = Math.min(H - 1, (y0i + R));
    const xStart = Math.max(0, (x0i - R));
    const xEnd   = Math.min(W - 1, (x0i + R));

    for (let y = yStart; y <= yEnd; y++) {
        const wyv = wy[(y - y0i) + R];
        const row = y * W;
        for (let x = xStart; x <= xEnd; x++) {
            sigA[row + x] += amp * (wx[(x - x0i) + R] * wyv);
        }
    }
}

// Read scanlines checkbox robustly from opts (UI uses id "cb-scanlines").
function getScanlinesFlag(opts) {
    if (!opts) return false;

    // common naming patterns
    const direct =
        (opts.cb_scanlines ?? opts.cbScanlines ?? opts.scanlines_flag ?? opts.scanlines ?? opts.enable_scanlines);
    if (typeof direct === 'boolean') return direct;

    // sometimes UI values are nested
    const ui = opts.ui || opts.controls || opts.flags || null;
    if (ui) {
        const v = (ui.cb_scanlines ?? ui.cbScanlines ?? ui.scanlines_flag ?? ui.scanlines ?? ui.enable_scanlines);
        if (typeof v === 'boolean') return v;
    }
    return false;
}

// ---------------------------
// Main: STM render
// ---------------------------


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

        // clipping
        low_clip = null,
        high_clip = null,

        // "slice" controls (education depth)
        focal_z = 0.0,
        hide_front = false,

        // scale bar (optional, same pattern as TEM)
        show_scale_bar = false,
        scale_bar_corner = 'bl',
        scale_bar_margin_px = 12,

        // canvas ctx (optional)
        canvasCtx = null
    } = opts;

    const [H, W] = img_size;

    // If empty, still clear canvas to background.
    if (!atoms || !atoms.length) {
        const out0 = new Uint8ClampedArray(H * W);
        out0.fill(clamp255(background_gray));
        if (canvasCtx) {
            const imageData = canvasCtx.createImageData(W, H);
            for (let i = 0, p = 0; i < out0.length; i++, p += 4) {
                const v = out0[i];
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
                    invert
                });
            }
        }
        return out0;
    }

    const [coords, scale, z_view] = compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel, camera);

    // Base STM "current map" signal (positive float image)
    // We'll build into sig.a (Float32), background 0.
    const sig = newFloatImage(H, W, 0);
    const sigA = sig.a;

    // Internal STM tuning (no UI changes required)
    // z falloff: larger alpha -> more surface-like
    const stm_alpha = (opts.stm_alpha != null) ? Number(opts.stm_alpha) : 2.8; // 1/Å
    const stm_sigma_mul = (opts.stm_sigma_mul != null) ? Number(opts.stm_sigma_mul) : 0.75;
    const stm_sigma_exp = (opts.stm_sigma_exp != null) ? Number(opts.stm_sigma_exp) : 1.05;

    // signal strength (maps normalized sig -> pixel delta)
    const stm_signal_strength = (opts.stm_signal_strength != null) ? Number(opts.stm_signal_strength) : 190.0;

    // scanline amps (in normalized signal units)
    const scan_off_amp = (opts.stm_scan_off_amp != null) ? Number(opts.stm_scan_off_amp) : 0.020;
    const scan_gain_amp = (opts.stm_scan_gain_amp != null) ? Number(opts.stm_scan_gain_amp) : 0.060;
    const scan_band_amp = (opts.stm_scan_band_amp != null) ? Number(opts.stm_scan_band_amp) : 0.016;
    const scan_drift_amp_px = (opts.stm_scan_drift_amp_px != null) ? Number(opts.stm_scan_drift_amp_px) : 1.10;

    // slow-variation smoothing (in rows)
    const scan_smooth_sigma = (opts.stm_scan_smooth_sigma != null) ? Number(opts.stm_scan_smooth_sigma) : 10.0;

    // Build signal: per-atom ROI Gaussians with z-dependent weights.
    // NOTE: bonds disabled by design.
    const Zs = atoms.map(a => a.Z | 0);
    const radiiA = Zs.map(Z => get_covalent_radius(Z));

    const typical_bond_A = opts.typical_bond_A ?? 1.40;
    const cap_rel_bond = opts.stm_sigma_cap_rel_bond ?? 0.45;
    const cap_abs_px = opts.stm_sigma_cap_abs_px ?? 18.0;

    for (let i = 0; i < atoms.length; i++) {
        const p = coords[i];
        if (!p) continue;

        const x0 = p[0];
        const y0 = p[1];

        const zAtom = (z_view ? (z_view[i] ?? 0) : (atoms[i].z ?? 0));
        if (hide_front && (zAtom > focal_z + 1e-9)) continue;

        // z-dependent weight (surface-like):
        //  - if hide_front: only atoms behind focal plane remain; deeper atoms contribute less
        //  - else: symmetric falloff around focal_z
        const dz = hide_front ? Math.max(0, (focal_z - zAtom)) : Math.abs(zAtom - focal_z);
        const wZ = Math.exp(-stm_alpha * dz);

        // 2D gaussian size (px): scaled covalent radius, but slightly smaller than TEM
        const rA = radiiA[i];

        // Per-element overrides (STM: apply size only)
        const sym = _EO_atomSymbol(atoms[i]);
        const [sizeMul, _darkMul] = _EO_getMul(ov, sym);
        let sigma = Math.max(1e-6, Math.pow(Math.max(1e-6, rA), stm_sigma_exp) * scale * stm_sigma_mul);
        sigma *= sizeMul;
        // keep a small minimum so distant zoom doesn't disappear completely
        if (sigma < 0.35) sigma = 0.35;

        // LOD cap to avoid huge blobs on heavy atoms / dense crystals
        const cap_px = Math.min(cap_abs_px, Math.max(1.0, cap_rel_bond * typical_bond_A * scale));
        if (sigma > cap_px) sigma = cap_px;

        // amplitude: weak Z dependence (STM is not pure Z, but gives variety)
        const Z = Math.max(1, Zs[i]);
        const aZ = Math.pow(Z, 0.35);

        const amp = wZ * aZ;
        if (amp < 1e-6) continue;

        addGaussianROI(sigA, H, W, x0, y0, sigma, amp);
    }

    // Compress dynamic range a bit (fast, stable)
    for (let p = 0; p < sigA.length; p++) sigA[p] = Math.sqrt(sigA[p]);

    // Normalize signal to ~[0..1]
    let maxSig = 0;
    for (let p = 0; p < sigA.length; p++) if (sigA[p] > maxSig) maxSig = sigA[p];
    const invMax = 1.0 / Math.max(1e-9, maxSig);
    for (let p = 0; p < sigA.length; p++) sigA[p] *= invMax;

    // Optional scanline block (true row artifacts + gentle drift)
    const scanlines_on = getScanlinesFlag(opts);
    if (scanlines_on) {
        // per-row offsets / gains / drift (low-frequency via smooth random walk)
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

        // low-frequency banding (slow wave over y)
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

            // drift (subpixel row shift)
            shiftRowLinear(sigA, rowOff, rowBuf, W, dx[y]);

            // apply gain + offset + banding
            const g = 1.0 + gain[y];
            const o = off[y] + band;

            for (let x = 0; x < W; x++) {
                let v = rowBuf[x] * g + o;
                if (v < 0) v = 0;
                sigA[rowOff + x] = v;
            }
        }

        // Keep within [0..1] again (after row mods)
        let max2 = 0;
        for (let p = 0; p < sigA.length; p++) if (sigA[p] > max2) max2 = sigA[p];
        const inv2 = 1.0 / Math.max(1e-9, max2);
        for (let p = 0; p < sigA.length; p++) sigA[p] *= inv2;
    }

    // Convert signal -> grayscale float image (img.a) with background reference.
    // STM is typically displayed as bumps/holes on a background; we map signal into darkening.
    const img = newFloatImage(H, W, background_gray);
    for (let p = 0; p < img.a.length; p++) {
        img.a[p] = background_gray - (stm_signal_strength * sigA[p]);
    }

    // Postprocess (matches TEM order)
    const scratch = new Float32Array(H * W);

    if (blur_sigma > 0) gaussianBlurFloat(img, blur_sigma, scratch);

    if (noise_stddev > 0) {
        const σ = noise_stddev;
        for (let i = 0; i < img.a.length; i++) img.a[i] += randn() * σ;
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

    // to 8-bit
    const out = new Uint8ClampedArray(H * W);
    for (let i = 0; i < out.length; i++) out[i] = clamp255(img.a[i]);

    // RENDER to canvas (grayscale -> RGBA)
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
                invert
            });
        }
    }

    return out;
}

// Backward-compatible alias (some code may call render_image_stm)
export function render_image_stm(atoms, opts = {}) {
    return render_stm_like(atoms, opts);
}
