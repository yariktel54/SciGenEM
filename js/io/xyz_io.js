// js/io/xyz_io.js
// Minimal XYZ loader (client-side). Bonds do NOT exist in XYZ => bonds=null.
//
// XYZ format:
//  line1: N
//  line2: comment/title
//  next N lines: element x y z  (Å)

import { fallback_symbol_to_Z, norm_sym } from '../chem/ptable_fallback.js';

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

function looksLikeRawXyzText(s) {
  s = String(s || '');
  if (!s.includes('\n')) return false;
  const lines = s.split(/\r?\n/);
  const n = parseInt(String(lines[0] || '').trim(), 10);
  if (!Number.isFinite(n) || n <= 0) return false;
  return lines.length >= (n + 2);
}

function guessTitleFromPath(pathOrText) {
  const s = String(pathOrText || '');
  if (!s) return 'XYZ';
  // handle blob urls with fragments: blob:...#file.xyz
  const frag = s.split('#').pop() || s;
  const base = frag.split('/').pop() || frag;
  return base || 'XYZ';
}

export async function load_xyz(pathOrText) {
  let text = pathOrText ?? '';
  const s = String(text);

  // If it's a URL/blob, fetch it unless it's clearly raw XYZ text.
  if (!looksLikeRawXyzText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
    const r = await fetch(s);
    if (!r.ok) throw new Error('XYZ: не вдалося завантажити (' + r.status + ')');
    text = await r.text();
  }

  const raw = String(text || '');
  const lines0 = raw.split(/\r?\n/);
  // trim empty edges
  while (lines0.length && String(lines0[0] || '').trim() === '') lines0.shift();
  while (lines0.length && String(lines0[lines0.length - 1] || '').trim() === '') lines0.pop();
  if (lines0.length < 3) throw new Error('XYZ: файл занадто короткий.');

  const n = parseInt(String(lines0[0] || '').trim(), 10);
  if (!Number.isFinite(n) || n <= 0) throw new Error('XYZ: перший рядок має бути числом N.');
  if (lines0.length < n + 2) throw new Error('XYZ: очікується ' + n + ' рядків атомів (отримано ' + Math.max(0, lines0.length - 2) + ').');

  const comment = String(lines0[1] || '').trim();
  const title = comment || guessTitleFromPath(pathOrText);

  const atoms = [];
  for (let i = 0; i < n; i++) {
    const line = String(lines0[2 + i] || '').trim();
    if (!line) throw new Error('XYZ: порожній рядок атома #' + (i + 1));
    const parts = line.split(/\s+/).filter(Boolean);
    if (parts.length < 4) throw new Error('XYZ: рядок має вигляд "El x y z": ' + line);
    const sym = parts[0];
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);
    if (![x, y, z].every(Number.isFinite)) throw new Error('XYZ: некоректні координати: ' + line);

    // Support numeric atomic number in XYZ (e.g. "6 0.0 0.0 0.0").
    let Z = 0;
    if (/^\d+$/.test(sym)) {
      Z = parseInt(sym, 10) | 0;
    } else {
      Z = symbol_to_Z(sym);
    }
    if (!Number.isFinite(Z) || Z <= 0) throw new Error('XYZ: невідомий елемент: ' + sym);

    atoms.push({ Z: Z | 0, x: +x, y: +y, z: +z });
  }

  return {
    title,
    kind: 'molecule',
    atoms,
    bonds: null
  };
}
