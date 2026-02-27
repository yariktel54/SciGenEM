// js/io/json_io.js
// Minimal JSON loader (client-side).
//
// Minimal schema:
// {
//   "title":"...",
//   "kind":"molecule"|"periodic",
//   "cell":{"a":[...],"b":[...],"c":[...]} (optional),
//   "atoms":[{"element":"C" or "Z":6, "x":..,"y":..,"z":..}, ...],
//   "bonds": null or [[i,j,type], ...] (optional)
// }
//
// Bonds semantics:
//  - missing "bonds" => bonds=null
//  - bonds=null => unknown (guess allowed later, not here)
//  - bonds=[] or array => explicit (guess forbidden)

import { fallback_symbol_to_Z, norm_sym } from '../chem/ptable_fallback.js';
import { RDKit } from '../chem/rdkit_wrap.js';
import { cell_vectors, frac_to_cart } from '../system/lattice.js';

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

function safe_symbol_to_Z(sym) {
  // Try RDKit periodic table first (same as CIF path), then fall back to local tables.
  const n = norm_sym(sym);
  let Z = 0;
  try { Z = RDKit.ptable_getAtomicNumber(n) || 0; } catch (_) { Z = 0; }
  if (!Z) Z = fallback_symbol_to_Z(n) || 0;
  if (!Z) Z = Z_BY_SYM[n] || 0;
  return Z | 0;
}


function looksLikeRawJsonText(s) {
  s = String(s || '').replace(/^\uFEFF/, '').trim();
  if (!s) return false;
  const c0 = s[0];
  if (!(c0 === '{' || c0 === '[')) return false;
  return true;
}

function guessTitleFromPath(pathOrText) {
  const s = String(pathOrText || '');
  if (!s) return 'JSON';
  const frag = s.split('#').pop() || s;
  const base = frag.split('/').pop() || frag;
  return base || 'JSON';
}

function dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function norm(a) { return Math.sqrt(dot(a,a)); }
function angle_deg(u, v) {
  const nu = norm(u), nv = norm(v);
  if (!(nu > 0 && nv > 0)) return 90;
  let c = dot(u, v) / (nu * nv);
  if (c > 1) c = 1;
  if (c < -1) c = -1;
  return Math.acos(c) * 180.0 / Math.PI;
}

function cell_vectors_to_params(cell) {
  if (!cell || typeof cell !== 'object') return null;
  const a = cell.a, b = cell.b, c = cell.c;
  if (!Array.isArray(a) || !Array.isArray(b) || !Array.isArray(c)) return null;
  if (a.length < 3 || b.length < 3 || c.length < 3) return null;
  const av = [Number(a[0]), Number(a[1]), Number(a[2])];
  const bv = [Number(b[0]), Number(b[1]), Number(b[2])];
  const cv = [Number(c[0]), Number(c[1]), Number(c[2])];
  if (![...av, ...bv, ...cv].every(Number.isFinite)) return null;
  const alen = norm(av);
  const blen = norm(bv);
  const clen = norm(cv);
  const alpha = angle_deg(bv, cv); // ∠(b,c)
  const beta  = angle_deg(av, cv); // ∠(a,c)
  const gamma = angle_deg(av, bv); // ∠(a,b)
  return { params: [alen, blen, clen, alpha, beta, gamma], vectors: { a: av, b: bv, c: cv } };
}

function cell_matrix_to_params(matrix) {
  if (!Array.isArray(matrix) || matrix.length < 3) return null;
  const a = matrix[0], b = matrix[1], c = matrix[2];
  if (!Array.isArray(a) || !Array.isArray(b) || !Array.isArray(c)) return null;
  if (a.length < 3 || b.length < 3 || c.length < 3) return null;
  const av = [Number(a[0]), Number(a[1]), Number(a[2])];
  const bv = [Number(b[0]), Number(b[1]), Number(b[2])];
  const cv = [Number(c[0]), Number(c[1]), Number(c[2])];
  if (![...av, ...bv, ...cv].every(Number.isFinite)) return null;
  const alen = norm(av);
  const blen = norm(bv);
  const clen = norm(cv);
  const alpha = angle_deg(bv, cv); // ∠(b,c)
  const beta  = angle_deg(av, cv); // ∠(a,c)
  const gamma = angle_deg(av, bv); // ∠(a,b)
  return { params: [alen, blen, clen, alpha, beta, gamma], vectors: { a: av, b: bv, c: cv } };
}

function best_species_element(speciesArr) {
  // Pymatgen Structure JSON: species is an array of {element, occu, ...}
  if (!Array.isArray(speciesArr) || speciesArr.length === 0) return null;
  let best = null;
  let bestOccu = -1;
  for (let i = 0; i < speciesArr.length; i++) {
    const s = speciesArr[i];
    if (!s) continue;
    const occu = (s.occu == null) ? 1 : Number(s.occu);
    const el = s.element || s.symbol || s.species || null;
    if (!el) continue;
    if (!Number.isFinite(occu)) {
      if (best == null) best = el;
      continue;
    }
    if (occu > bestOccu) {
      bestOccu = occu;
      best = el;
    }
  }
  return best;
}

function parse_pymatgen_structure(obj, titleFallback) {
  // Recognize pymatgen Structure.as_dict() JSON.
  // Keys: "@class":"Structure", lattice:{matrix/a/b/c/alpha/beta/gamma,pbc}, sites:[{species:[...], abc:[...], xyz:[...]}]
  if (!obj || typeof obj !== 'object') return null;
  const cls = obj['@class'] || obj['@CLASS'] || obj.class || null;
  if (String(cls || '') !== 'Structure') return null;
  if (!obj.lattice || (!obj.lattice.matrix && (obj.lattice.a == null)) || !Array.isArray(obj.sites)) return null;

  const title = (obj.title != null && String(obj.title).trim()) ? String(obj.title).trim() : titleFallback;

  const pbc = obj.lattice.pbc;
  const kind = (Array.isArray(pbc) && pbc.some(Boolean)) ? 'periodic' : 'molecule';

  // --- Cell params (CIF-style): [a,b,c,alpha,beta,gamma] ---
  let a = Number(obj.lattice.a), b = Number(obj.lattice.b), c = Number(obj.lattice.c);
  let alpha = Number(obj.lattice.alpha), beta = Number(obj.lattice.beta), gamma = Number(obj.lattice.gamma);

  if (![a, b, c, alpha, beta, gamma].every(Number.isFinite)) {
    const cm = cell_matrix_to_params(obj.lattice.matrix);
    if (!cm) throw new Error('JSON: lattice.matrix некоректний.');
    a = cm.params[0]; b = cm.params[1]; c = cm.params[2];
    alpha = cm.params[3]; beta = cm.params[4]; gamma = cm.params[5];
  }

  const cellParams = [a, b, c, alpha, beta, gamma];

  // IMPORTANT: to match CIF pipeline, we intentionally build the conventional basis from (a,b,c,α,β,γ),
  // and convert fractional "abc" -> cart using this basis.
  const vv = cell_vectors(a, b, c, alpha, beta, gamma);
  const vx = vv[0], vy = vv[1], vz = vv[2];
  const cellVectors = { a: vx, b: vy, c: vz };

  const atoms_frac = [];
  const atoms_cart = [];

  for (let i = 0; i < obj.sites.length; i++) {
    const site = obj.sites[i];
    if (!site || typeof site !== 'object') throw new Error("JSON: sites[" + i + "] не є об'єктом.");

    const sym = best_species_element(site.species) || site.element || site.symbol || null;
    if (!sym) throw new Error("JSON: sites[" + i + "] не має елемента (species).");
    const Z = safe_symbol_to_Z(sym);
    if (!Number.isFinite(Z) || Z <= 0) throw new Error("JSON: невідомий елемент sites[" + i + "] = " + sym);

    const label = (site.label != null) ? String(site.label) : '';

    // Prefer fractional coords for periodic structures (pymatgen always provides 'abc').
    if (Array.isArray(site.abc) && site.abc.length >= 3) {
      const fx = Number(site.abc[0]);
      const fy = Number(site.abc[1]);
      const fz = Number(site.abc[2]);
      if (![fx, fy, fz].every(Number.isFinite)) throw new Error("JSON: sites[" + i + "].abc некоректні.");

      atoms_frac.push({ Z: Z | 0, fx: +fx, fy: +fy, fz: +fz, label });

      const p = frac_to_cart([fx, fy, fz], vx, vy, vz);
      atoms_cart.push({ Z: Z | 0, x: p[0], y: p[1], z: p[2], label });
      continue;
    }

    // Fallback: cartesian xyz (mostly for non-periodic dumps)
    if (Array.isArray(site.xyz) && site.xyz.length >= 3) {
      const x = Number(site.xyz[0]);
      const y = Number(site.xyz[1]);
      const z = Number(site.xyz[2]);
      if (![x, y, z].every(Number.isFinite)) throw new Error("JSON: sites[" + i + "].xyz некоректні.");
      atoms_cart.push({ Z: Z | 0, x: +x, y: +y, z: +z, label });
      continue;
    }

    throw new Error("JSON: sites[" + i + "] не має abc/xyz координат.");
  }

  // Pymatgen dict doesn't include bonds.
  return {
    title,
    kind,
    cell: cellParams,
    cell_from: 'json',
    cell_vectors: cellVectors,
    atoms: atoms_cart,
    atoms_cart,
    atoms_frac,
    bonds: null
  };
}



function parse_bond_order(t) {
  if (t == null) return 1;
  if (typeof t === 'number') {
    const o = (t | 0);
    if ([1,2,3,4].includes(o)) return o;
    return 1;
  }
  const s = String(t).trim().toLowerCase();
  if (!s) return 1;
  if (s === '1' || s === 'single' || s === 's') return 1;
  if (s === '2' || s === 'double' || s === 'd') return 2;
  if (s === '3' || s === 'triple' || s === 't') return 3;
  if (s === '4' || s === 'aromatic' || s === 'a') return 4;
  const n = parseInt(s, 10);
  if ([1,2,3,4].includes(n)) return n;
  return 1;
}

export async function load_json_system(pathOrText) {
  let text = pathOrText ?? '';
  const s = String(text);

  // If it's a URL/blob, fetch it unless it's clearly raw JSON text.
  if (!looksLikeRawJsonText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
    const r = await fetch(s);
    if (!r.ok) throw new Error('JSON: не вдалося завантажити (' + r.status + ')');
    text = await r.text();
  }

  const raw = String(text || '').replace(/^\uFEFF/, '').trim();
  if (!raw) throw new Error('JSON: порожній вміст.');

  let obj;
  try {
    obj = JSON.parse(raw);
  } catch (e) {
    throw new Error('JSON: помилка парсингу.');
  }

  // Allow "atoms" provided as the root array.
  if (Array.isArray(obj)) obj = { atoms: obj };
  if (!obj || typeof obj !== 'object') throw new Error('JSON: очікується об\'єкт або масив.');

  const titleFallback = guessTitleFromPath(pathOrText);

  // Support pymatgen Structure JSON out-of-the-box.
  const pmg = parse_pymatgen_structure(obj, titleFallback);
  if (pmg) return pmg;

  const title = (obj.title != null && String(obj.title).trim()) ? String(obj.title).trim() : titleFallback;
  const kind = (obj.kind === 'periodic' || obj.kind === 'molecule') ? obj.kind : 'molecule';

  // atoms
  const atomsIn = obj.atoms;
  if (!Array.isArray(atomsIn) || atomsIn.length === 0) throw new Error('JSON: atoms має бути непорожнім масивом.');

  const atoms = [];
  for (let i = 0; i < atomsIn.length; i++) {
    const a = atomsIn[i];
    if (!a || typeof a !== 'object') throw new Error('JSON: atoms[' + i + '] не є об\'єктом.');

    let Z = 0;
    if ('Z' in a) {
      Z = parseInt(a.Z, 10);
      if (!Number.isFinite(Z) || Z <= 0) throw new Error('JSON: atoms[' + i + '].Z некоректний.');
    } else if ('element' in a) {
      Z = safe_symbol_to_Z(a.element);
      if (!Number.isFinite(Z) || Z <= 0) throw new Error('JSON: невідомий елемент atoms[' + i + '].element=' + a.element);
    } else {
      throw new Error('JSON: atoms[' + i + '] має містити element або Z.');
    }

    const x = parseFloat(a.x);
    const y = parseFloat(a.y);
    const z = ('z' in a) ? parseFloat(a.z) : 0.0;
    if (![x, y, z].every(Number.isFinite)) throw new Error('JSON: atoms[' + i + '] координати некоректні.');

    atoms.push({ Z: Z | 0, x: +x, y: +y, z: +z });
  }

  // bonds
  let bonds;
  if (!('bonds' in obj)) {
    bonds = null;
  } else if (obj.bonds === null) {
    bonds = null;
  } else if (Array.isArray(obj.bonds)) {
    bonds = [];
    const N = atoms.length;
    for (let k = 0; k < obj.bonds.length; k++) {
      const b = obj.bonds[k];
      if (!Array.isArray(b) || b.length < 2) throw new Error('JSON: bonds[' + k + '] має бути [i,j,type].');
      const i = parseInt(b[0], 10);
      const j = parseInt(b[1], 10);
      if (![i, j].every(Number.isFinite)) throw new Error('JSON: bonds[' + k + '] індекси некоректні.');
      if (!(i >= 0 && i < N && j >= 0 && j < N) || i === j) throw new Error('JSON: bonds[' + k + '] індекси поза межами або i==j.');
      const o = parse_bond_order(b[2]);
      bonds.push([i, j, o]);
    }
  } else {
    throw new Error('JSON: bonds має бути null або масивом.');
  }

  // cell (optional)
  let cellParams = null;
  let cellVectors = null;
  if (obj.cell) {
    const cv = cell_vectors_to_params(obj.cell);
    if (cv) {
      cellParams = cv.params;
      cellVectors = cv.vectors;
    }
  }

  const atoms_cart = atoms;
  const atoms_frac = [];

  return {
    title,
    kind,
    cell: cellParams,
    cell_from: 'json',
    cell_vectors: cellVectors,
    atoms: atoms_cart,
    atoms_cart,
    atoms_frac,
    bonds
  };
}