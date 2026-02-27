// js/io/poscar_io.js
// Minimal POSCAR/CONTCAR (VASP5) loader (client-side). Bonds do NOT exist => bonds=null.
//
// Supported layout (VASP5):
//  1: comment/title
//  2: scale (float; if <0 => target volume)
//  3-5: 3 lattice vectors (3 floats each)
//  6: element symbols (e.g. "U As Se")
//  7: counts (e.g. "1 1 1")
//  8: optional "Selective dynamics" line
//  8/9: "Direct" or "Cartesian" (first letter D/C)
//  then N coordinate lines (N = sum(counts)); extra selective flags are ignored.

import { fallback_symbol_to_Z, norm_sym } from '../chem/ptable_fallback.js';
import { frac_to_cart } from '../system/lattice.js';

// Full symbol→Z fallback (covers elements not included in ptable_fallback.js)
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

function symbol_to_Z(sym) {
  const n = norm_sym(sym);
  let z = fallback_symbol_to_Z(n);
  if (z) return z;
  z = Z_BY_SYM[n] || 0;
  return z;
}

function dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function cross(a, b) {
  return [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  ];
}
function norm(v) { return Math.sqrt(Math.max(0, dot(v, v))); }

function angle_deg(u, v) {
  const nu = norm(u);
  const nv = norm(v);
  if (nu < 1e-12 || nv < 1e-12) return 90;
  let c = dot(u, v) / (nu * nv);
  if (c > 1) c = 1;
  if (c < -1) c = -1;
  return (Math.acos(c) * 180.0) / Math.PI;
}

function vectors_to_cell(v1, v2, v3) {
  const a = norm(v1);
  const b = norm(v2);
  const c = norm(v3);
  const alpha = angle_deg(v2, v3);
  const beta = angle_deg(v1, v3);
  const gamma = angle_deg(v1, v2);
  return [a, b, c, alpha, beta, gamma];
}

function build_cart_to_frac(vx, vy, vz) {
  const bxc = cross(vy, vz);
  const cxa = cross(vz, vx);
  const axb = cross(vx, vy);
  const det = dot(vx, bxc);
  if (!Number.isFinite(det) || Math.abs(det) < 1e-12) return null;
  const inv = [
    [bxc[0] / det, bxc[1] / det, bxc[2] / det],
    [cxa[0] / det, cxa[1] / det, cxa[2] / det],
    [axb[0] / det, axb[1] / det, axb[2] / det]
  ];
  return function cart_to_frac(v) {
    return [
      inv[0][0]*v[0] + inv[0][1]*v[1] + inv[0][2]*v[2],
      inv[1][0]*v[0] + inv[1][1]*v[1] + inv[1][2]*v[2],
      inv[2][0]*v[0] + inv[2][1]*v[1] + inv[2][2]*v[2]
    ];
  };
}

function looksLikeRawPoscarText(s) {
  s = String(s || '');
  if (!s.includes('\n')) return false;
  const lines = s.split(/\r?\n/).map(x => String(x || '').trim()).filter(x => x !== '');
  if (lines.length < 8) return false;
  const scale = parseFloat(lines[1]);
  if (!Number.isFinite(scale) || scale === 0) return false;
  // 3 lattice vector lines: 3 numbers each
  for (let i = 2; i <= 4; i++) {
    const p = lines[i].split(/\s+/).filter(Boolean);
    if (p.length < 3) return false;
    const v = [parseFloat(p[0]), parseFloat(p[1]), parseFloat(p[2])];
    if (!v.every(Number.isFinite)) return false;
  }
  // VASP5: symbols line then counts
  const symT = lines[5].split(/\s+/).filter(Boolean);
  const cntT = lines[6].split(/\s+/).filter(Boolean);
  if (!symT.length || !cntT.length) return false;
  // require at least one symbol that looks like element
  if (!symT.every(t => /^[A-Za-z]{1,2}$/.test(t))) return false;
  if (!cntT.every(t => /^-?\d+$/.test(t))) return false;
  return true;
}

function guessTitleFromPath(pathOrText, fallbackTitleLine) {
  const t = String(fallbackTitleLine || '').trim();
  if (t) return t;
  const s = String(pathOrText || '');
  if (!s) return 'POSCAR';
  const frag = s.split('#').pop() || s;
  const base = frag.split('/').pop() || frag;
  return base || 'POSCAR';
}

export async function load_poscar_system(pathOrText) {
  let text = pathOrText ?? '';
  const s = String(text);

  // If it's a URL/blob, fetch it unless it is clearly raw POSCAR text.
  if (!looksLikeRawPoscarText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
    const r = await fetch(s);
    if (!r.ok) throw new Error('POSCAR: не вдалося завантажити (' + r.status + ')');
    text = await r.text();
  }

  const raw = String(text || '');
  const lines0 = raw.split(/\r?\n/);
  // trim empty edges
  while (lines0.length && String(lines0[0] || '').trim() === '') lines0.shift();
  while (lines0.length && String(lines0[lines0.length - 1] || '').trim() === '') lines0.pop();
  if (lines0.length < 8) throw new Error('POSCAR: файл занадто короткий.');

  const titleLine = String(lines0[0] || '').trim();
  const title = guessTitleFromPath(pathOrText, titleLine);

  const scaleRaw = parseFloat(String(lines0[1] || '').trim());
  if (!Number.isFinite(scaleRaw) || scaleRaw === 0) throw new Error('POSCAR: некоректний scale у рядку 2.');

  function parse_vec(line, idx) {
    const p = String(line || '').trim().split(/\s+/).filter(Boolean);
    if (p.length < 3) throw new Error('POSCAR: некоректний вектор ґратки (рядок ' + idx + ').');
    const v = [parseFloat(p[0]), parseFloat(p[1]), parseFloat(p[2])];
    if (!v.every(Number.isFinite)) throw new Error('POSCAR: некоректні числа у векторі ґратки (рядок ' + idx + ').');
    return v;
  }

  let v1 = parse_vec(lines0[2], 3);
  let v2 = parse_vec(lines0[3], 4);
  let v3 = parse_vec(lines0[4], 5);

  // Apply scale factor to lattice vectors (and Cartesian coords).
  // If scaleRaw < 0: interpret as target volume (Å^3).
  let scale = scaleRaw;
  if (scaleRaw < 0) {
    const targetV = Math.abs(scaleRaw);
    const curV = Math.abs(dot(v1, cross(v2, v3)));
    if (curV < 1e-12) throw new Error('POSCAR: нульовий об’єм ґратки.');
    scale = Math.cbrt(targetV / curV);
  }

  v1 = [v1[0] * scale, v1[1] * scale, v1[2] * scale];
  v2 = [v2[0] * scale, v2[1] * scale, v2[2] * scale];
  v3 = [v3[0] * scale, v3[1] * scale, v3[2] * scale];

  const cell = vectors_to_cell(v1, v2, v3);
  const vx = v1, vy = v2, vz = v3;
  const cart_to_frac = build_cart_to_frac(vx, vy, vz);

  const symTokens = String(lines0[5] || '').trim().split(/\s+/).filter(Boolean);
  const cntTokens = String(lines0[6] || '').trim().split(/\s+/).filter(Boolean);
  if (!symTokens.length || !cntTokens.length) throw new Error('POSCAR: очікуються рядки елементів і кількостей (VASP5).');

  // VASP5 requirement: symbols line must not be purely numeric.
  if (symTokens.every(t => /^-?\d+$/.test(t))) {
    throw new Error('POSCAR: формат без рядка символів (VASP4) не підтримується (потрібен VASP5).');
  }

  const symbols = symTokens.map(norm_sym);
  const counts = cntTokens.map(x => parseInt(x, 10) | 0);
  if (symbols.length !== counts.length) {
    throw new Error('POSCAR: кількість символів елементів не збігається з кількістю чисел у рядку counts.');
  }
  if (!counts.every(n => Number.isFinite(n) && n >= 0)) {
    throw new Error('POSCAR: некоректні counts.');
  }

  const nAtoms = counts.reduce((a, b) => a + (b | 0), 0);
  if (nAtoms <= 0) throw new Error('POSCAR: nAtoms=0.');

  let iLine = 7;
  const l7 = String(lines0[iLine] || '').trim();
  if (/^s/i.test(l7)) {
    // Selective dynamics
    iLine++;
  }

  const modeLine = String(lines0[iLine] || '').trim();
  if (!modeLine) throw new Error('POSCAR: відсутній рядок Direct/Cartesian.');
  const c0 = modeLine[0].toLowerCase();
  const isDirect = (c0 === 'd');
  const isCart = (c0 === 'c' || c0 === 'k');
  if (!isDirect && !isCart) throw new Error('POSCAR: очікується Direct або Cartesian.');
  iLine++;

  const atoms_cart = [];
  const atoms_frac = [];

  let specieIndex = 0;
  let specieLeft = counts[0] | 0;
  let Zcur = symbol_to_Z(symbols[0]);
  if (!Zcur) throw new Error('POSCAR: невідомий елемент: ' + symbols[0]);

  function stepSpecies() {
    while (specieIndex < counts.length && specieLeft <= 0) {
      specieIndex++;
      if (specieIndex >= counts.length) break;
      specieLeft = counts[specieIndex] | 0;
      Zcur = symbol_to_Z(symbols[specieIndex]);
      if (!Zcur) throw new Error('POSCAR: невідомий елемент: ' + symbols[specieIndex]);
    }
  }

  for (let i = 0; i < nAtoms; i++) {
    stepSpecies();
    if (specieIndex >= counts.length) throw new Error('POSCAR: забагато координатних рядків.');

    const line = String(lines0[iLine + i] || '').trim();
    if (!line) throw new Error('POSCAR: порожній рядок координат #' + (i + 1));
    const p = line.split(/\s+/).filter(Boolean);
    if (p.length < 3) throw new Error('POSCAR: координати мають містити щонайменше 3 числа.');
    const a0 = parseFloat(p[0]);
    const a1 = parseFloat(p[1]);
    const a2 = parseFloat(p[2]);
    if (![a0, a1, a2].every(Number.isFinite)) throw new Error('POSCAR: некоректні координати: ' + line);

    if (isDirect) {
      const fx = +a0, fy = +a1, fz = +a2;
      atoms_frac.push({ Z: Zcur | 0, fx, fy, fz });
      const xyz = frac_to_cart([fx, fy, fz], vx, vy, vz);
      atoms_cart.push({ Z: Zcur | 0, x: xyz[0], y: xyz[1], z: xyz[2] });
    } else {
      // Cartesian coordinates are scaled by the same factor as lattice vectors.
      const x = +a0 * scale;
      const y = +a1 * scale;
      const z = +a2 * scale;
      atoms_cart.push({ Z: Zcur | 0, x, y, z });
      if (typeof cart_to_frac === 'function') {
        const f = cart_to_frac([x, y, z]);
        atoms_frac.push({ Z: Zcur | 0, fx: f[0], fy: f[1], fz: f[2] });
      }
    }

    specieLeft--;
  }

  return {
    title,
    kind: 'periodic',
    cell,
    cell_from: 'poscar',
    atoms_cart,
    atoms_frac: atoms_frac.length ? atoms_frac : null,
    // convenience (contract expects cart)
    atoms: atoms_cart,
    bonds: null
  };
}
