// js/render/render_tem.js
// (split from renderer.js, Wave R1)

import { get_covalent_radius, guess_bonds_by_distance } from '../bonds/bonds.js';
import {
    compute_scaled_coordinates,
    newFloatImage,
    clamp255,
    gaussianBlurFloat,
    getGaussian1D_cached,
    draw_gaussian_line,
    draw_bond_glow_line,
    draw_scale_bar
} from './render_common.js';


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

export function draw_atoms(img, atoms, coords, scale, background_gray, opts = {}) {
    const {
        compose = 'sum',
        focal_z = 0.0,
        dof_strength = 0.0,
        hide_front = false
    } = opts;

    // (compose/dof_strength поки не використовуються напряму — залишено для сумісності)
    void compose; void dof_strength;

    const Zs = atoms.map(a => a.Z | 0);
    const radiiA = Zs.map(Z => get_covalent_radius(Z));

    // Element overrides (per-element multipliers)
    const ov = (opts && opts.element_overrides) ? opts.element_overrides : null;

    // v3: 2D Gaussian як wx * wy (1D ваги по X/Y)
    // v4: wx/wy беруться з кешу getGaussian1D_cached() (менше exp() + майже без алокацій)

    for (let i = 0; i < atoms.length; i++) {
        const x0 = coords[i][0];
        const y0 = coords[i][1];

        const zv = (opts && Array.isArray(opts.z_view)) ? opts.z_view : null;
        if (hide_front && ((zv ? (zv[i] ?? 0) : (atoms[i].z ?? 0)) > focal_z + 1e-9)) continue;

        const sym = _EO_atomSymbol(atoms[i]);
        const [sizeMul, darkMul] = _EO_getMul(ov, sym);

        // Реалістична "геометрична" ширина атома: залежить від ковалентного радіуса
        // atom_size_mul — загальний множник, atom_size_exp — підсилює різницю між H і важкими
        const atom_size_mul = opts.atom_size_mul ?? 1.2;     // 0.6..1.6 (тюнінг)
        const atom_size_exp = opts.atom_size_exp ?? 1.25;    // 1.0..1.5 (1.25 дає помітну різницю)
        const rA = radiiA[i];

        // без "мінімумів": лише числовий захист від нуля
        // + LOD-кап по масштабу, щоб важкі елементи/кристали не перетворювались на "чорну пляму"
        //   (кап прив’язаний до типового bond-length у пікселях, а не до формату)
        let sigma = Math.max(1e-6, Math.pow(rA, atom_size_exp) * scale * atom_size_mul);
        sigma *= sizeMul;
        const typical_bond_A = opts.typical_bond_A ?? 1.40;               // Å
        const cap_rel_bond = opts.atom_sigma_cap_rel_bond ?? 0.55;        // частка від bondPx
        const cap_abs_px = opts.atom_sigma_cap_abs_px ?? 24.0;            // абсолютний safety-cap
        const cap_px = Math.min(cap_abs_px, Math.max(1.5, cap_rel_bond * typical_bond_A * scale));
        if (sigma > cap_px) sigma = cap_px;

        // Інтенсивність
        const atom_dark_mul = opts.atom_dark_mul ?? 2.5;   // 1.5..4.0
        const atom_dark_exp = opts.atom_dark_exp ?? 1.2;   // 1.1..1.6
        const baseZ = Math.pow(Zs[i], atom_dark_exp) * darkMul;
        const intensity = -(baseZ * atom_dark_mul);

        // «ядро»
        const core_sigma_rel = opts.core_sigma_rel ?? 0.18; // 0.12..0.22
        const core_rel = opts.core_rel ?? 0.6;             // 0.3..0.9
        const core_Z0 = opts.core_Z0 ?? 12;                // 8..20

        const core_sigma = sigma * core_sigma_rel;

        // приглушення ядра для великих Z
        const core_falloff = 1.0 / (1.0 + (Zs[i] / core_Z0) * (Zs[i] / core_Z0));
        const core_intensity = baseZ * core_rel * core_falloff;

        // Опорна ціла точка (subpixel): t = x - floor(x0), frac = x0 - floor(x0)
        // Це дозволяє рахувати 2D Gaussian як wx[t]*wy[u] без exp() у внутрішньому циклі.
        const x0i = Math.floor(x0);
        const y0i = Math.floor(y0);
        const fracX = x0 - x0i;
        const fracY = y0 - y0i;

        // --- Gaussian "тіло" атома (subpixel) ---
        const R = Math.ceil(3 * sigma);

        // v2: якщо атом повністю за кадром — пропускаємо весь атом
        if ((x0 + R) < 0 || (x0 - R) >= img.w || (y0 + R) < 0 || (y0 - R) >= img.h) continue;

        // v3: 1D Gaussian ваги (subpixel)
        const wx = getGaussian1D_cached(sigma, fracX, R);
        const wy = getGaussian1D_cached(sigma, fracY, R);

        // v2: обрізаємо межі обходу одразу (менше if у циклі)
        const yStart = Math.max(0, (y0i - R));
        const yEnd   = Math.min(img.h - 1, (y0i + R));
        const xStart = Math.max(0, (x0i - R));
        const xEnd   = Math.min(img.w - 1, (x0i + R));

        for (let y = yStart; y <= yEnd; y++) {
            const wyv = wy[(y - y0i) + R];
            for (let x = xStart; x <= xEnd; x++) {
                const v = intensity * (wx[(x - x0i) + R] * wyv);
                const p = (y * img.w + x);
                const current = img.a[p];
                const candidate = background_gray + v;
                img.a[p] = Math.min(current, candidate);
            }
        }

        // --- «ядро» (яскравий центр) (subpixel) ---
        const Rc = Math.ceil(3 * core_sigma);

        // v2: ядро може бути поза кадром навіть коли тіло частково видно
        if (!((x0 + Rc) < 0 || (x0 - Rc) >= img.w || (y0 + Rc) < 0 || (y0 - Rc) >= img.h)) {
            // v3: окремі 1D ваги для ядра (інший sigma)
            const wxC = getGaussian1D_cached(core_sigma, fracX, Rc);
            const wyC = getGaussian1D_cached(core_sigma, fracY, Rc);

            const yStartC = Math.max(0, (y0i - Rc));
            const yEndC   = Math.min(img.h - 1, (y0i + Rc));
            const xStartC = Math.max(0, (x0i - Rc));
            const xEndC   = Math.min(img.w - 1, (x0i + Rc));

            for (let y = yStartC; y <= yEndC; y++) {
                const wyv = wyC[(y - y0i) + Rc];
                for (let x = xStartC; x <= xEndC; x++) {
                    const add = core_intensity * (wxC[(x - x0i) + Rc] * wyv);
                    const p = (y * img.w + x);
                    img.a[p] += add;
                }
            }
        }
    }
}

// opts: {wave_width_px, wave_amp}
export function draw_bonds(img, coords, bonds, atoms, opts = {}) {
    const { wave_width_px = 8.0, wave_amp = 1.0 } = opts;
    const halfw = Math.max(1.0, wave_width_px * 0.5);

    for (const b of bonds) {
        const [i, j, bt] = b;
        const p1 = coords[i], p2 = coords[j];
        if (!p1 || !p2) continue;

        const Z1 = atoms[i]?.Z ?? 0;
        const Z2 = atoms[j]?.Z ?? 0;
        const Zavg = 0.5 * (Z1 + Z2);

        const base_int = Math.min(255.0, Math.pow(Zavg, 0.9));
        const glow_int = Math.min(255.0, Math.pow(Zavg, 1.2) * 2.0);
        const bond_intensity = Math.min(255.0, base_int * wave_amp);
        const glow_intensity = Math.min(255.0, glow_int * wave_amp);

        // зсуви для 1x/2x/3x
        const direction = [(p2[0] - p1[0]), (p2[1] - p1[1])];
        const L = Math.hypot(direction[0], direction[1]) || 1;
        const nrm = [-direction[1] / L, direction[0] / L];
        const offset_base = Math.max(1.0, 0.8 * halfw);

        // v2 bbox-кулінг: якщо сегмент + товщина повністю поза кадром — пропускаємо
        const max_off = (bt === 2) ? offset_base : (bt === 3) ? (1.2 * offset_base) : 0.0;
        const pad = halfw + 2 + max_off;
        const minx = Math.min(p1[0], p2[0]) - pad;
        const maxx = Math.max(p1[0], p2[0]) + pad;
        const miny = Math.min(p1[1], p2[1]) - pad;
        const maxy = Math.max(p1[1], p2[1]) + pad;
        if (maxx < 0 || minx >= img.w || maxy < 0 || miny >= img.h) continue;


        let offsets = [[0, 0]];
        if (bt === 2) offsets = [[nrm[0] * offset_base, nrm[1] * offset_base], [-nrm[0] * offset_base, -nrm[1] * offset_base]];
        else if (bt === 3) offsets = [[nrm[0] * 1.2 * offset_base, nrm[1] * 1.2 * offset_base], [0, 0], [-nrm[0] * 1.2 * offset_base, -nrm[1] * 1.2 * offset_base]];

        for (const off of offsets) {
            const q1 = [p1[0] + off[0], p1[1] + off[1]];
            const q2 = [p2[0] + off[0], p2[1] + off[1]];
            draw_gaussian_line(img, q1, q2, bond_intensity, halfw);

            // glow тільки для центральної лінії
            if (off[0] === 0 && off[1] === 0) {
                draw_bond_glow_line(img, p1, p2, glow_intensity * 0.85, Math.max(0.5, 0.6 * halfw));
            }
        }
    }
}

// ---------------------------
// AFM-like surface render (Wave 2)
//  - Height/topography map (max of Gaussian bumps)
//  - Edge enhancement (simple Laplacian)
//  - Optional single-line "bond contrast" overlay (bt ignored; legacy profile with glow; no multi-lines)
//  - NOTE: DoF is intentionally ignored in AFM mode (surface mode).
// ---------------------------


export function render_image_tem(
    atoms,
    {
        bonds = null,
        img_size = [400, 400],
        angstroms_per_pixel = 0.1,
        blur_sigma = 2,
        background_gray = 127,
        invert = false,
        noise_stddev = 2.0,
        contrast = 1.0,
        compose_mode = 'sum',
        draw_bonds_flag = true,
        // camera/viewState for stable centering + pan/rotate
        camera = null,
        bond_wave_width_px = 8.0,
        bond_wave_amplitude = 1.0,
        low_clip = null,
        high_clip = null,
        focal_z = 0.0,
        dof_strength = 0.0,
        hide_front = false,
        show_scale_bar = false,
        element_overrides = null,
        scale_bar_corner = 'bl',
        scale_bar_margin_px = 12,

        // --- v5: тюнінг атомів (опційно; якщо не передано — використовуються дефолти draw_atoms) ---
        atom_size_mul = 1.2,
        atom_size_exp = 1.25,
        atom_dark_mul = 2.5,
        atom_dark_exp = 1.2,
        core_sigma_rel = 0.18,
        core_rel = 0.6,
        core_Z0 = 12,

        // --- v5: FOV prep (опційно) ---
        // Якщо true — перед DoF/зв'язками робимо легкий кулінг атомів/зв’язків, які гарантовано поза кадром.
        // Це НЕ змінює геометрію, тільки зменшує обсяг роботи, коли структура велика/зум сильний.
        fov_prep = true,

        // додатково: контекст Canvas для малювання масштабної лінійки поверх
        canvasCtx = null
    } = {}
) {
    const [H, W] = img_size;
    const img = newFloatImage(H, W, background_gray);

    if (!atoms || !atoms.length) {
        // навіть якщо пусто — очистимо канву
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

    const [coords, scale, z_view] = compute_scaled_coordinates(atoms, img_size, angstroms_per_pixel, camera);

    // Семантика bonds:
    //  • bonds === null або undefined → зв'язки невідомі → можна здогадувати по відстанях;
    //  • bonds — масив (навіть порожній [] ) → це явна топологія, НІЧОГО не "guess"-ити.
    let use_bonds = [];
    if (draw_bonds_flag) {
        if (bonds == null) {
            // bonds unknown -> guess allowed, but we MUST avoid heavy O(N^2) on large/zoomed-out views.
            const MAX_GUESS_N = 6000;
            const typicalBondA = 1.4; // rough typical covalent bond length
            const bondPxEst = typicalBondA / Math.max(angstroms_per_pixel, 1e-9);
            const minVisiblePx = Math.max(3.0, bond_wave_width_px * 0.35);
            if (atoms.length <= MAX_GUESS_N && bondPxEst >= minVisiblePx) {
                use_bonds = guess_bonds_by_distance(atoms);
                try { use_bonds.__tem_guessed = true; } catch (_) {}
            } else {
                use_bonds = []; // skip guessing for this view
            }
        } else {
            use_bonds = bonds;
        }
    }

    // LOD для bonds: якщо в цьому масштабі зв’язки майже не видно — не малюємо їх
    // ВАЖЛИВО: LOD застосовуємо лише до guessed-bonds. Для explicit bonds (MOL/SMILES/CSV з BONDS)
    // користувач очікує, що вони будуть видимі при увімкненому чекбоксі.
    const bondsAreGuessed = (bonds == null) || (!!(bonds && bonds.__tem_guessed));
    let drawBondsEffective = (draw_bonds_flag && use_bonds && use_bonds.length > 0);
    if (drawBondsEffective && bondsAreGuessed) {
        const maxSamples = 256;
        const step = Math.max(1, Math.floor(use_bonds.length / maxSamples));
        let sum = 0, cnt = 0;

        for (let k = 0; k < use_bonds.length; k += step) {
            const b = use_bonds[k];
            const i = b[0], j = b[1];
            const p1 = coords[i], p2 = coords[j];
            if (!p1 || !p2) continue;
            const dx = p2[0] - p1[0];
            const dy = p2[1] - p1[1];
            const L = Math.hypot(dx, dy);
            if (Number.isFinite(L)) { sum += L; cnt++; }
            if (cnt >= maxSamples) break;
        }

        const avgLenPx = cnt ? (sum / cnt) : 0;
        const minVisiblePx = Math.max(3.0, bond_wave_width_px * 0.35);
        if (avgLenPx < minVisiblePx) drawBondsEffective = false;
    }
    if (!drawBondsEffective) use_bonds = [];

    // --- v5: FOV pre-pass (атомний) ---
    // Потрібно насамперед для DoF-гілки (щоб не робити K×N по всіх атомах/зв'язках).
    let atoms_draw = atoms;
    let coords_draw = coords;
    let z_draw = z_view;
    let bonds_draw = use_bonds;

    if (fov_prep) {
        const Zs = atoms.map(a => a.Z | 0);
        const radiiA = Zs.map(Z => get_covalent_radius(Z));

        const idxs = [];
        idxs.length = 0;

        for (let i = 0; i < atoms.length; i++) {
            const x0 = coords[i][0];
            const y0 = coords[i][1];

            // hide_front: відразу викидаємо «передні» атоми з prep (еквівалентно draw_atoms)
            if (hide_front && ((z_view ? (z_view[i] ?? 0) : (atoms[i].z ?? 0)) > focal_z + 1e-9)) continue;

            const rA = radiiA[i];
            let sigma = Math.max(1e-6, Math.pow(rA, atom_size_exp) * scale * atom_size_mul);
            // те ж саме LOD-обмеження, що й у draw_atoms()
            const typical_bond_A = 1.40;
            const cap_rel_bond = 0.55;
            const cap_abs_px = 24.0;
            const cap_px = Math.min(cap_abs_px, Math.max(1.5, cap_rel_bond * typical_bond_A * scale));
            if (sigma > cap_px) sigma = cap_px;
            const R = Math.ceil(3 * sigma);

            // повністю поза кадром -> не малюємо
            if ((x0 + R) < 0 || (x0 - R) >= W || (y0 + R) < 0 || (y0 - R) >= H) continue;

            idxs.push(i);
        }

        // якщо майже все в кадрі — не робимо зайвих копій
        if (idxs.length > 0 && idxs.length < atoms.length * 0.98) {
            const map = new Int32Array(atoms.length);
            map.fill(-1);

            atoms_draw = new Array(idxs.length);
            coords_draw = new Array(idxs.length);
            z_draw = z_view ? new Array(idxs.length) : null;

            for (let li = 0; li < idxs.length; li++) {
                const gi = idxs[li];
                map[gi] = li;
                atoms_draw[li] = atoms[gi];
                coords_draw[li] = coords[gi];
                if (z_draw) z_draw[li] = z_view[gi];
            }

            if (drawBondsEffective && bonds_draw && bonds_draw.length) {
                const filtered = [];
                for (let k = 0; k < bonds_draw.length; k++) {
                    const b = bonds_draw[k];
                    const gi = b[0], gj = b[1];
                    const li = map[gi], lj = map[gj];
                    if (li < 0 || lj < 0) continue;
                    filtered.push([li, lj, b[2]]);
                }
                bonds_draw = filtered;
            } else {
                bonds_draw = [];
            }
        } else {
            // нічого не змінюємо
            bonds_draw = drawBondsEffective ? bonds_draw : [];
        }
    } else {
        bonds_draw = drawBondsEffective ? bonds_draw : [];
    }

    const atomDrawOpts = {
        compose: compose_mode,
        focal_z,
        dof_strength: 0,
        hide_front,
        z_view: z_draw,
        element_overrides,

        // тюнінги атомів
        atom_size_mul,
        atom_size_exp,
        atom_dark_mul,
        atom_dark_exp,
        core_sigma_rel,
        core_rel,
        core_Z0
    };

    const dof_on = (dof_strength && dof_strength > 1e-6);

    function getZDraw(i) {
        return z_draw ? (z_draw[i] ?? 0) : (atoms_draw[i].z ?? 0);
    }

    // scratch буфер для GaussianBlur (уникнути великих алокацій)
    const scratch = new Float32Array(H * W);

    if (dof_on) {
        // --- v5: формуємо шари за 1 прохід ---
        // zmin/zmax беремо з *atoms_draw* (бо позакадрове не має значення)
        let zmin = Infinity, zmax = -Infinity;
        for (let i = 0; i < atoms_draw.length; i++) {
            const z = getZDraw(i);
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        if (!Number.isFinite(zmin) || !Number.isFinite(zmax)) { zmin = 0; zmax = 0; }

        const K = (zmax > zmin) ? 12 : 1;
        const inv = (zmax > zmin) ? (K / (zmax - zmin)) : 0;

        const layers = Array.from({ length: K }, () => []);
        for (let i = 0; i < atoms_draw.length; i++) {
            const z = getZDraw(i);
            let bi = 0;
            if (K > 1) {
                bi = Math.floor((z - zmin) * inv);
                if (bi < 0) bi = 0;
                if (bi >= K) bi = K - 1;
            }
            layers[bi].push(i);
        }

        // bonds: НЕ рендеримо всередині DoF-шарів (інакше DoF/blur «з’їдає» тонкі зв’язки).
        // Bonds будуть намальовані окремим оверлеєм ПІСЛЯ DoF-композиції.

        // Буфери шару (reuse)
        const layer = { a: new Float32Array(H * W), h: H, w: W, bg: background_gray };
        const tmp = { a: new Float32Array(H * W), h: H, w: W, bg: 0 };

        for (let bi = 0; bi < K; bi++) {
            const z_c = (K === 1) ? (0.5 * (zmin + zmax)) : (zmin + (bi + 0.5) * (zmax - zmin) / K);
            if (hide_front && z_c > focal_z + 1e-9) continue;

            const idxs = layers[bi];
            const bsub = null; // bonds rendered post-DoF overlay

            if ((!idxs || !idxs.length)) continue;

            // локальні масиви (для draw_atoms потрібна щільна індексація)
            const atoms_sub = new Array(idxs.length);
            const coords_sub = new Array(idxs.length);
            const indexMap = new Int32Array(atoms_draw.length);
            indexMap.fill(-1);

            for (let li = 0; li < idxs.length; li++) {
                const gi = idxs[li];
                indexMap[gi] = li;
                atoms_sub[li] = atoms_draw[gi];
                coords_sub[li] = coords_draw[gi];
            }

            // render layer
            layer.a.fill(background_gray);
            draw_atoms(layer, atoms_sub, coords_sub, scale, background_gray, { ...atomDrawOpts, focal_z: 0, hide_front: false });
// tmp = layer - bg
            for (let p = 0; p < tmp.a.length; p++) tmp.a[p] = layer.a[p] - background_gray;

            // DoF blur sigma
            const sigma_px = Math.abs(z_c - focal_z) * dof_strength / Math.max(angstroms_per_pixel, 1e-9);
            if (sigma_px > 0.35) gaussianBlurFloat(tmp, sigma_px, scratch);

            // accumulate
            for (let p = 0; p < tmp.a.length; p++) img.a[p] += tmp.a[p];
        }
    } else {
        // no DoF
        draw_atoms(img, atoms_draw, coords_draw, scale, background_gray, atomDrawOpts);
        if (drawBondsEffective && bonds_draw && bonds_draw.length) {
            draw_bonds(img, coords_draw, bonds_draw, atoms_draw, { wave_width_px: bond_wave_width_px, wave_amp: bond_wave_amplitude });
        }
    }

    if (blur_sigma > 0) gaussianBlurFloat(img, blur_sigma, scratch);

    if (dof_on && drawBondsEffective && bonds_draw && bonds_draw.length) {
        // --- Bonds overlay (Option A) ---
        // Малюємо зв’язки ПІСЛЯ DoF-композиції, щоб вони не «розчинялись» у DoF blur.
        // На bonds застосовуємо 0 або дуже малий blur (не більше 0.6 px),
        // щоб вони лишались помітними при blur_sigma=1..3.
        const bondLayer = { a: new Float32Array(H * W), h: H, w: W, bg: 0 };
        draw_bonds(bondLayer, coords_draw, bonds_draw, atoms_draw, {
            wave_width_px: bond_wave_width_px,
            wave_amp: bond_wave_amplitude
        });

        // Раніше blur для bonds був занадто малий, і зв’язки виглядали «надто різкими».
        // Тут робимо м’який blur ближчий до загального, але з верхнім капом,
        // щоб bonds не «розчинялись» при великих blur_sigma.
        const bonds_blur = (blur_sigma > 0)
            ? Math.min(1.4, Math.max(0.0, blur_sigma) * 0.65)
            : 0.0;
        if (bonds_blur > 0.35) gaussianBlurFloat(bondLayer, bonds_blur, scratch);

        for (let p = 0; p < bondLayer.a.length; p++) img.a[p] += bondLayer.a[p];
    }

    if (noise_stddev > 0) {
        const σ = noise_stddev;
        for (let i = 0; i < img.a.length; i++) {
            // Box–Muller
            const u1 = Math.random(), u2 = Math.random();
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
        const a = Math.min(lo, hi), b = Math.max(lo, hi);
        for (let i = 0; i < img.a.length; i++) {
            const v = img.a[i];
            img.a[i] = v < a ? a : v > b ? b : v;
        }
    }

    // до 8-біт
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
