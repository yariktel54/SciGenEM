// js/io/cif_io.js
// Мінімальний парсер CIF як у Python, але з фолбеками.
// load_cif: повертає { cell, atoms_frac, atoms_cart, bonds }
//
// Зміни відносно «старого» файлу:
//  - НЕ кидаємо помилку, якщо немає _cell_length_a/b/c (часто в молекулярних CIF).
//  - Підтримуємо і fractional (_atom_site_fract_*), і cartesian (_atom_site_Cartn_* / _atom_site_cartn_*).
//  - Підтримуємо _geom_bond_* як опційні bonds (інакше bonds=null).

import { RDKit } from '../chem/rdkit_wrap.js';
import { norm_sym, fallback_symbol_to_Z } from '../chem/ptable_fallback.js';

// Full symbol→Z fallback (covers heavy elements even when RDKit WASM isn't ready)
const ALL_SYMBOLS = [
  null,
  'H','He',
  'Li','Be','B','C','N','O','F','Ne',
  'Na','Mg','Al','Si','P','S','Cl','Ar',
  'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
  'Ga','Ge','As','Se','Br','Kr',
  'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
  'In','Sn','Sb','Te','I','Xe',
  'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
  'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
  'Tl','Pb','Bi','Po','At','Rn',
  'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
  'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn',
  'Nh','Fl','Mc','Lv','Ts','Og'
];
const Z_BY_SYM = (() => {
  const m = Object.create(null);
  for (let i = 1; i < ALL_SYMBOLS.length; i++) m[ALL_SYMBOLS[i]] = i;
  return m;
})();


function safe_symbol_to_Z(sym) {
    const n = norm_sym(sym);
    let Z = 0;
    try { Z = RDKit.ptable_getAtomicNumber(n) || 0; } catch { Z = 0; }
    if (!Z) Z = fallback_symbol_to_Z(n) || 0;
    if (!Z) Z = Z_BY_SYM[n] || 0;
    return Z || 6;
}

function looksLikeRawCifText(s) {
    s = String(s || '');
    if (!s.includes('\n')) return false;
    const t = s.trimStart();
    return t.startsWith('data_') || t.startsWith('loop_') || t.startsWith('_');
}

function grabNumToken(str) {
    // 12.34(5) -> 12.34; "12.3" -> 12.3
    let s = String(str ?? '').trim();
    s = s.replace(/^['"]|['"]$/g, '');
    s = s.replace(/\([^\)]*\)/g, '');
    const m = s.match(/[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?/);
    if (!m) return NaN;
    return parseFloat(m[0]);
}

function bboxCellFromCart(atoms_cart) {
    if (!atoms_cart || atoms_cart.length < 2) return null;
    let minx = Infinity, miny = Infinity, minz = Infinity;
    let maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
    for (const a of atoms_cart) {
        if (a.x < minx) minx = a.x;
        if (a.y < miny) miny = a.y;
        if (a.z < minz) minz = a.z;
        if (a.x > maxx) maxx = a.x;
        if (a.y > maxy) maxy = a.y;
        if (a.z > maxz) maxz = a.z;
    }
    const dx = Math.max(1e-6, maxx - minx);
    const dy = Math.max(1e-6, maxy - miny);
    const dz = Math.max(1e-6, maxz - minz);
    return [dx * 1.05, dy * 1.05, dz * 1.05, 90.0, 90.0, 90.0];
}

export async function load_cif(pathOrText) {
    // У браузері «path» = або URL/blob:, або прямий CIF-текст
    let text = pathOrText ?? "";

    const s = String(text);
    if (!looksLikeRawCifText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
        const r = await fetch(s);
        if (!r.ok) throw new Error("CIF: не вдалося завантажити (" + r.status + ")");
        text = await r.text();
    }

    const lines = String(text)
        .split(/\r?\n/)
        .map(s => s.trim())
        .filter(Boolean);

    if (!lines.length) throw new Error("CIF: порожній вхід");

    function grabFloat(key) {
        key = key.toLowerCase();
        for (const ln of lines) {
            if (ln.toLowerCase().startsWith(key)) {
                const last = ln.split(/\s+/).at(-1);
                const v = grabNumToken(last);
                if (!Number.isNaN(v)) return v;
            }
        }
        return null;
    }

    // Cell (може бути відсутній)
    const a = grabFloat('_cell_length_a');
    const b = grabFloat('_cell_length_b');
    const c = grabFloat('_cell_length_c');
    const alpha = grabFloat('_cell_angle_alpha') ?? 90.0;
    const beta  = grabFloat('_cell_angle_beta') ?? 90.0;
    const gamma = grabFloat('_cell_angle_gamma') ?? 90.0;

    const atoms_frac = [];
    const atoms_cart = [];

    // Для bonds по label
    const labelToIndex_frac = new Map();
    const labelToIndex_cart = new Map();
    const pendingBondLabels = []; // [label1,label2]

    for (let i = 0; i < lines.length; i++) {
        if (lines[i].toLowerCase().startsWith('loop_')) {
            i++;
            const headers = [];
            while (i < lines.length && lines[i].startsWith('_')) { headers.push(lines[i]); i++; }
            const lower = headers.map(h => h.toLowerCase());

            const col = (names) => {
                for (const cand of names) {
                    const k = lower.indexOf(cand.toLowerCase());
                    if (k >= 0) return k;
                }
                return -1;
            };

            // ---------- ATOMS LOOP ----------
            if (lower.some(h => h.startsWith('_atom_site_'))) {
                const ix = col(['_atom_site_fract_x']);
                const iy = col(['_atom_site_fract_y']);
                const iz = col(['_atom_site_fract_z']);

                const jx = col(['_atom_site_cartn_x', '_atom_site_Cartn_x']);
                const jy = col(['_atom_site_cartn_y', '_atom_site_Cartn_y']);
                const jz = col(['_atom_site_cartn_z', '_atom_site_Cartn_z']);

                const ilabel = col(['_atom_site_label']);
                const isym = col(['_atom_site_type_symbol', '_atom_site_label']);

                const hasFrac = (ix >= 0 && iy >= 0 && iz >= 0);
                const hasCart = (jx >= 0 && jy >= 0 && jz >= 0);

                while (i < lines.length && !lines[i].startsWith('_') && !lines[i].toLowerCase().startsWith('loop_')) {
                    const parts = lines[i].split(/\s+/);
                    i++;
                    const need = Math.max(ix, iy, iz, jx, jy, jz, ilabel, isym);
                    if (need < 0 || parts.length <= need) continue;

                    const symRaw = parts[isym].replace(/^['"]|['"]$/g, '');
                    const m = /^([A-Za-z]+)/.exec(symRaw);
                    const sym = m ? m[1] : symRaw;
                    const Z = safe_symbol_to_Z(sym);

                    const label = (ilabel >= 0 ? parts[ilabel] : '').replace(/^['"]|['"]$/g, '').trim();

                    if (hasFrac) {
                        const fx = grabNumToken(parts[ix]);
                        const fy = grabNumToken(parts[iy]);
                        const fz = grabNumToken(parts[iz]);
                        if (![fx, fy, fz].some(Number.isNaN)) {
                            const idx = atoms_frac.length;
                            atoms_frac.push({ Z, fx, fy, fz, label });
                            if (label) labelToIndex_frac.set(label, idx);
                        }
                    }

                    if (hasCart) {
                        const x = grabNumToken(parts[jx]);
                        const y = grabNumToken(parts[jy]);
                        const z = grabNumToken(parts[jz]);
                        if (![x, y, z].some(Number.isNaN)) {
                            const idx = atoms_cart.length;
                            atoms_cart.push({ Z, x, y, z, label });
                            if (label) labelToIndex_cart.set(label, idx);
                        }
                    }
                }

                i--; // компенсуємо outer for
                continue;
            }

            // ---------- BONDS LOOP (optional) ----------
            if (lower.some(h => h.startsWith('_geom_bond_'))) {
                const b1 = col(['_geom_bond_atom_site_label_1']);
                const b2 = col(['_geom_bond_atom_site_label_2']);
                if (b1 >= 0 && b2 >= 0) {
                    while (i < lines.length && !lines[i].startsWith('_') && !lines[i].toLowerCase().startsWith('loop_')) {
                        const parts = lines[i].split(/\s+/);
                        i++;
                        const need = Math.max(b1, b2);
                        if (parts.length <= need) continue;
                        pendingBondLabels.push([parts[b1], parts[b2]]);
                    }
                    i--; // компенсуємо
                }
                continue;
            }
        }
    }

    if (!atoms_frac.length && !atoms_cart.length) {
        throw new Error("CIF: не знайдено _atom_site_fract_* або _atom_site_Cartn_*");
    }

    // cell decision + provenance
    //  - 'cif'  : true periodic metadata present
    //  - 'bbox' : fallback for molecular CIF with Cartn coords
    //  - 'unit' : fallback unit cube for fractional-only molecular CIF
    let cell = null;
    let cell_from = 'none';
    if (a != null && b != null && c != null) {
        cell = [a, b, c, alpha, beta, gamma];
        cell_from = 'cif';
    } else if (atoms_cart.length > 0) {
        cell = bboxCellFromCart(atoms_cart);
        cell_from = 'bbox';
    } else if (atoms_frac.length > 0) {
        cell = [1, 1, 1, 90.0, 90.0, 90.0];
        cell_from = 'unit';
    }

    // bonds mapping by label (optional)
    let bonds = null;
    if (pendingBondLabels.length) {
        const map = (atoms_frac.length > 0) ? labelToIndex_frac : labelToIndex_cart;
        const out = [];
        for (const [r1, r2] of pendingBondLabels) {
            const l1 = String(r1 || '').replace(/^['"]|['"]$/g, '').trim();
            const l2 = String(r2 || '').replace(/^['"]|['"]$/g, '').trim();
            const iA = map.get(l1);
            const iB = map.get(l2);
            if (iA == null || iB == null || iA === iB) continue;
            out.push([iA, iB, 1]);
        }
        bonds = out.length ? out : null;
    }

    return { cell, cell_from, atoms_frac, atoms_cart, bonds };
}
