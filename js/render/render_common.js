// js/render/render_common.js
// (split from renderer.js, Wave R1)

// ---------------------------
// V4: кеш 1D Gaussian ваг
// ---------------------------
// У v3 exp() перенесено з внутрішнього піксельного циклу у побудову 1D-ваг.
// Тут робимо наступний крок: 1D ваги кешуються за (R, sigma, subpixel-фаза),
// щоб для сотень/тисяч атомів одного типу не будувати wx/wy заново.
// Для максимальної візуальної стабільності sigma квантуємо дуже тонко (1e-3 px).

const _GAUSS1D_CACHE = new Map();
const _GAUSS1D_CACHE_MAX = 4096; // LRU по insertion order Map
const _SIGMA_Q = 1000; // 1e-3 px
const _FRAC_BINS = 64; // subpixel фази

function _gaussKey(R, sigmaBin, fracBin) {
    return R + '|' + sigmaBin + '|' + fracBin;
}

export function getGaussian1D_cached(sigma, frac, R) {
    // sigmaBin: 0.001 px; fracBin: 1/64
    const sigmaBin = Math.max(1, Math.round(sigma * _SIGMA_Q));
    let fracBin = Math.round(frac * _FRAC_BINS);
    if (fracBin < 0) fracBin = 0;
    // важливо: frac ∈ [0,1). При round() інколи може вилізти рівно _FRAC_BINS,
    // тоді краще притиснути до останнього біну, а не робити fracQ=1.0.
    if (fracBin >= _FRAC_BINS) fracBin = _FRAC_BINS - 1;

    const key = _gaussKey(R, sigmaBin, fracBin);
    let w = _GAUSS1D_CACHE.get(key);
    if (w) {
        // refresh LRU: delete+set
        _GAUSS1D_CACHE.delete(key);
        _GAUSS1D_CACHE.set(key, w);
        return w;
    }

    const sigmaQ = sigmaBin / _SIGMA_Q;
    const fracQ = fracBin / _FRAC_BINS;
    const inv2sig2 = 1.0 / (2.0 * sigmaQ * sigmaQ);
    w = new Float32Array(2 * R + 1);
    for (let t = -R; t <= R; t++) {
        const d = t - fracQ;
        w[t + R] = Math.exp(-(d * d) * inv2sig2);
    }

    _GAUSS1D_CACHE.set(key, w);
    if (_GAUSS1D_CACHE.size > _GAUSS1D_CACHE_MAX) {
        const firstKey = _GAUSS1D_CACHE.keys().next().value;
        if (firstKey !== undefined) _GAUSS1D_CACHE.delete(firstKey);
    }
    return w;
}

// camera (optional): { pan_px:[dx,dy], rotZ_rad, center_A:[cx,cy,cz] }
export function compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel, camera) {
    const scale = 1.0 / Math.max(angstroms_per_pixel, 1e-9);

    // стабільний центр кадру (для парних/непарних розмірів)
    const cx = (img_size[1] - 1) * 0.5;
    const cy = (img_size[0] - 1) * 0.5;

    // Центр моделі:
    //  - якщо є камера з center_A -> беремо її (стабільний anchor для lazy view)
    //  - інакше -> центр bounding-box (краще за mean для "A.B")
    let centerX = null;
    let centerY = null;
    if (camera && Array.isArray(camera.center_A)) {
        const v0 = Number(camera.center_A[0]);
        const v1 = Number(camera.center_A[1]);
        if (Number.isFinite(v0) && Number.isFinite(v1)) {
            centerX = v0;
            centerY = v1;
        }
    }

    if (centerX == null || centerY == null) {
        let xmin = Infinity, xmax = -Infinity;
        let ymin = Infinity, ymax = -Infinity;
        for (const a of atoms) {
            const x = a.x, y = a.y;
            if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }
        if (!Number.isFinite(xmin) || !Number.isFinite(xmax)) {
            xmin = xmax = 0;
            ymin = ymax = 0;
        }
        centerX = 0.5 * (xmin + xmax);
        centerY = 0.5 * (ymin + ymax);
    }

    const pan = (camera && Array.isArray(camera.pan_px)) ? camera.pan_px : [0, 0];
    const panX = Number(pan[0]) || 0;
    const panY = Number(pan[1]) || 0;

    // Rotations:
    //  - rotX/rotY: tilt model in 3D (mixes XY with Z)
    //  - rotZ: in-plane rotation (existing stage control)
    let angZ = 0;
    let angX = 0;
    let angY = 0;
    if (camera) {
        if (Number.isFinite(camera.rotZ_rad)) angZ = camera.rotZ_rad;
        else if (Number.isFinite(camera.rot_deg)) angZ = (camera.rot_deg * Math.PI) / 180.0;

        if (Number.isFinite(camera.rotX_rad)) angX = camera.rotX_rad;
        else if (Number.isFinite(camera.tilt_x_deg)) angX = (camera.tilt_x_deg * Math.PI) / 180.0;

        if (Number.isFinite(camera.rotY_rad)) angY = camera.rotY_rad;
        else if (Number.isFinite(camera.tilt_y_deg)) angY = (camera.tilt_y_deg * Math.PI) / 180.0;
    }

    const ca = Math.cos(angZ);
    const sa = Math.sin(angZ);
    const cxX = Math.cos(angX);
    const sxX = Math.sin(angX);
    const cyY = Math.cos(angY);
    const syY = Math.sin(angY);

    // centerZ (needed for tilt). Prefer camera.center_A[2], else bbox-center z.
    let centerZ = null;
    if (camera && Array.isArray(camera.center_A)) {
        const v2 = Number(camera.center_A[2]);
        if (Number.isFinite(v2)) centerZ = v2;
    }
    if (centerZ == null) {
        let zmin = Infinity, zmax = -Infinity;
        for (const a of atoms) {
            const z = (a.z ?? 0);
            if (!Number.isFinite(z)) continue;
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        if (!Number.isFinite(zmin) || !Number.isFinite(zmax)) {
            zmin = zmax = 0;
        }
        centerZ = 0.5 * (zmin + zmax);
    }

    const useTilt = (Math.abs(angX) > 1e-12) || (Math.abs(angY) > 1e-12);
    const z_view = useTilt ? new Array(atoms.length) : null;

    // ВАЖЛИВО: НЕ округляємо — subpixel координати
    const coords = new Array(atoms.length);
    for (let i = 0; i < atoms.length; i++) {
        const a = atoms[i];
        let x = (a.x - centerX);
        let y = (a.y - centerY);
        let z = ((a.z ?? 0) - centerZ);

        // tilt around X
        if (useTilt && Math.abs(angX) > 1e-12) {
            const y1 = y * cxX - z * sxX;
            const z1 = y * sxX + z * cxX;
            y = y1;
            z = z1;
        }
        // tilt around Y
        if (useTilt && Math.abs(angY) > 1e-12) {
            const x2 = x * cyY + z * syY;
            const z2 = -x * syY + z * cyY;
            x = x2;
            z = z2;
        }

        // in-plane rotation around Z
        const xr = x * ca - y * sa;
        const yr = x * sa + y * ca;

        coords[i] = [xr * scale + cx + panX, yr * scale + cy + panY];
        if (z_view) z_view[i] = z + centerZ;
    }

    return [coords, scale, z_view];
}



// «полотно» — масив Float32 (H×W), далі конвертуємо у 8-біт для Canvas
export function newFloatImage(h, w, bg = 127) {
    const a = new Float32Array(h * w);
    a.fill(bg);
    return { a, h, w, bg };
}
function idx(img, x, y) { return y * img.w + x; }
export function clamp255(v) { return v < 0 ? 0 : v > 255 ? 255 : v; }

function _stmSmoothstep(edge0, edge1, x) {
    if (!(edge1 > edge0)) return x >= edge1 ? 1 : 0;
    let t = (x - edge0) / (edge1 - edge0);
    if (t < 0) t = 0;
    else if (t > 1) t = 1;
    return t * t * (3 - 2 * t);
}

export function stm_orange_map_rgb(grayValue) {
    let t = clamp255(Number(grayValue)) / 255.0;
    t = Math.pow(t, 0.92);

    const orange = [255, 136, 28];
    const white = [255, 247, 236];

    const mid = _stmSmoothstep(0.06, 0.72, t);
    const hi = _stmSmoothstep(0.62, 1.00, t);

    let r = orange[0] * mid;
    let g = orange[1] * mid;
    let b = orange[2] * mid;

    r += (white[0] - r) * hi;
    g += (white[1] - g) * hi;
    b += (white[2] - b) * hi;

    return [clamp255(r) | 0, clamp255(g) | 0, clamp255(b) | 0];
}

export function draw_stm_orange_frame_to_canvas(canvasCtx, gray, w, h) {
    if (!canvasCtx) return;
    const imageData = canvasCtx.createImageData(w, h);
    const data = imageData.data;
    for (let i = 0, p = 0; i < gray.length; i++, p += 4) {
        const rgb = stm_orange_map_rgb(gray[i]);
        data[p] = rgb[0];
        data[p + 1] = rgb[1];
        data[p + 2] = rgb[2];
        data[p + 3] = 255;
    }
    canvasCtx.putImageData(imageData, 0, 0);
}

export function gaussianBlurFloat(img, sigma, scratch = null) {
    if (sigma <= 0) return img;
    const { h, w, a } = img;

    // scratch: Float32Array(h*w), використовується як проміжний буфер, щоб уникати великих алокацій
    let tmp = scratch;
    if (!(tmp instanceof Float32Array) || tmp.length !== a.length) {
        tmp = new Float32Array(a.length);
    }

    // kernel
    const ksize = Math.max(1, (Math.floor(sigma * 3) * 2 + 1));
    const half = ksize >> 1;
    const kernel = new Float32Array(ksize);
    let sum = 0;
    const inv2sig2 = 1.0 / (2 * sigma * sigma);
    for (let i = -half; i <= half; i++) {
        const v = Math.exp(-(i * i) * inv2sig2);
        kernel[i + half] = v;
        sum += v;
    }
    for (let i = 0; i < ksize; i++) kernel[i] /= sum;

    // X: a -> tmp
    for (let y = 0; y < h; y++) {
        const row = y * w;
        for (let x = 0; x < w; x++) {
            let s = 0;
            for (let k = -half; k <= half; k++) {
                const xx = (x + k < 0) ? 0 : (x + k >= w) ? (w - 1) : (x + k);
                s += a[row + xx] * kernel[k + half];
            }
            tmp[row + x] = s;
        }
    }

    // Y: tmp -> a
    for (let y = 0; y < h; y++) {
        const row = y * w;
        for (let x = 0; x < w; x++) {
            let s = 0;
            for (let k = -half; k <= half; k++) {
                const yy = (y + k < 0) ? 0 : (y + k >= h) ? (h - 1) : (y + k);
                s += tmp[yy * w + x] * kernel[k + half];
            }
            a[row + x] = s;
        }
    }

    return img;
}

// opts: {compose, focal_z, dof_strength, hide_front}

export function draw_gaussian_line(img, p1, p2, intensity = 30, halfwidth_px = 4.0) {
    const x1 = Math.round(p1[0]);
    const y1 = Math.round(p1[1]);
    const x2 = Math.round(p2[0]);
    const y2 = Math.round(p2[1]);
    let dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy) | 0;
    if (L <= 0) return;

    dx /= L; dy /= L;
    const nx = -dy, ny = dx;
    const radius = Math.ceil(halfwidth_px);

    for (let i = 0; i <= L; i++) {
        const x = Math.round(x1 + dx * i), y = Math.round(y1 + dy * i);
        for (let sx = -radius; sx <= radius; sx++) {
            for (let sy = -radius; sy <= radius; sy++) {
                const xx = x + sx, yy = y + sy;
                if (xx < 0 || xx >= img.w || yy < 0 || yy >= img.h) continue;
                const d_perp = Math.abs(sx * nx + sy * ny);
                const s = d_perp / Math.max(halfwidth_px, 1e-6);
                if (s <= 1.0) {
                    const w = Math.cos(0.5 * Math.PI * s);
                    const val = intensity * (w * w);
                    if (val > 0.25) img.a[idx(img, xx, yy)] += val;
                }
            }
        }
    }
}

export function draw_bond_glow_line(img, p1, p2, intensity = 50, sigma = 1.2) {
    const x1 = Math.round(p1[0]);
    const y1 = Math.round(p1[1]);
    const x2 = Math.round(p2[0]);
    const y2 = Math.round(p2[1]);
    let dx = x2 - x1, dy = y2 - y1;
    const L = Math.hypot(dx, dy) | 0;
    if (L === 0) return;

    for (let i = 0; i <= L; i++) {
        const x = Math.round(x1 + (dx * i) / L), y = Math.round(y1 + (dy * i) / L);
        const R = Math.max(1, Math.floor(3 * sigma));
        for (let sx = -R; sx <= R; sx++) {
            for (let sy = -R; sy <= R; sy++) {
                const xx = x + sx, yy = y + sy;
                if (xx < 0 || xx >= img.w || yy < 0 || yy >= img.h) continue;
                const dist2 = sx * sx + sy * sy;
                const value = intensity * Math.exp(-dist2 / (2 * sigma * sigma));
                if (value > 0.5) img.a[idx(img, xx, yy)] += value;
            }
        }
    }
}

function nice_scale_length(angstroms_per_pixel, width_px, min_px = 80, max_px = 200) {
    if (width_px <= 0 || angstroms_per_pixel <= 0) return [0, 0];
    const target_px = Math.min(Math.max(width_px * 0.25, min_px), max_px);
    const target_A = target_px * angstroms_per_pixel;
    const k = target_A > 0 ? Math.floor(Math.log10(target_A)) : 0;

    let best = null;
    for (let shift = -3; shift <= 3; shift++) {
        const base = Math.pow(10, k + shift);
        for (const m of [1, 2, 5]) {
            const valA = m * base;
            const valPx = valA / angstroms_per_pixel;
            if (valPx >= min_px && valPx <= max_px) {
                if (!best || Math.abs(valPx - target_px) < Math.abs(best[0] - target_px)) best = [valPx, valA];
            }
        }
    }
    if (!best) {
        const valPx = max_px, valA = valPx * angstroms_per_pixel;
        return [Math.round(valPx), valA];
    }
    return [Math.round(best[0]), best[1]];
}

export function draw_scale_bar(ctx, img, angstroms_per_pixel, { corner = 'bl', margin = 12, invert = false } = {}) {
    const h = img.h, w = img.w;
    const [bar_px, bar_A] = nice_scale_length(angstroms_per_pixel, w, 80, 200);
    if (bar_px <= 0) return;

    const col = invert ? 0 : 255;

    const thickness = Math.max(4, Math.min(18, (w / 80) | 0));
    const tick_w = Math.max(2, thickness - 2);
    const tick_h = (thickness * 1.8) | 0;
    const gap_txt = Math.max(8, (thickness >> 1) + 6);
    const pad_edge = Math.max(20, (margin * 1.6) | 0);

    const text = (Math.abs(bar_A - Math.round(bar_A)) < 1e-6)
        ? `${Math.round(bar_A)} A`
        : `${bar_A.toFixed(1)} A`;

    let x0, x1, y;
    const right = /r$/.test(corner);
    const top = /^t/.test(corner);

    if (right) {
        x1 = w - margin;
        x0 = x1 - bar_px;
    } else {
        x0 = margin;
        x1 = margin + bar_px;
    }

    y = top ? pad_edge : (h - pad_edge);

    ctx.fillStyle = `rgb(${col},${col},${col})`;

    ctx.fillRect(x0, y - (thickness >> 1), bar_px, thickness);
    ctx.fillRect(x0, y - (tick_h >> 1), tick_w, tick_h);
    ctx.fillRect(x1 - tick_w, y - (tick_h >> 1), tick_w, tick_h);

    ctx.font = `${Math.max(14, (thickness * 2) | 0)}px system-ui, Arial, sans-serif`;
    ctx.textBaseline = 'top';

    const tw = ctx.measureText(text).width | 0;
    const tx = right ? (x1 - tw) : x0;
    const ty = top
        ? (y + (thickness >> 1) + gap_txt)
        : (y - (thickness >> 1) - gap_txt - (thickness * 2));

    ctx.fillText(text, tx, Math.max(4, Math.min(h - 16, ty)));
}

