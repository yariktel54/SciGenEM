// js/bonds/bonds_guess.js
// Single bond-perception fallback (distance-based).
//
// Strict semantics (enforced by callers, NOT here):
//  - bonds === null  -> unknown -> guessing is allowed as a late fallback.
//  - bonds is array  -> explicit topology (even []) -> guessing is forbidden.
//
// This module provides ONE implementation to reuse across all formats.
//
// Notes:
//  - Atoms are expected in Å.
//  - For periodic structures, prefer running this on atoms that already include tiling/PBC
//    (i.e., atomsView from PeriodicProvider). The 'periodic/cell' opts are accepted for API
//    compatibility but are not required for the current pipeline.

import { fallback_symbol_to_Z, fallback_rcovalent, norm_sym } from '../chem/ptable_fallback.js';

// Compact full periodic table mapping for element->Z fallback
const _ALL_SYMBOLS = [
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

const _Z_BY_SYM = (() => {
  const m = Object.create(null);
  for (let i = 1; i < _ALL_SYMBOLS.length; i++) m[_ALL_SYMBOLS[i]] = i;
  return m;
})();

function _symbol_to_Z(sym) {
  const n = norm_sym(sym);
  return fallback_symbol_to_Z(n) || _Z_BY_SYM[n] || 0;
}

// Extra radii for common non-fallback elements (Å). Only used if fallback_rcovalent(Z) is 0.
const _RCOV_EXTRA = {
  21: 1.44, 22: 1.36, 23: 1.25, 24: 1.27, 25: 1.39,
  31: 1.22, 32: 1.20, 33: 1.19, 34: 1.20, 36: 1.16,
  37: 2.20, 38: 1.95, 39: 1.90, 40: 1.75, 41: 1.64, 42: 1.54,
  48: 1.44, 49: 1.42, 50: 1.39, 51: 1.39, 52: 1.38,
  54: 1.40,
  55: 2.44, 56: 2.15,
  57: 2.07, 58: 2.04, 59: 2.03, 60: 2.01, 61: 1.99, 62: 1.98, 63: 1.98, 64: 1.96,
  65: 1.94, 66: 1.92, 67: 1.92, 68: 1.89, 69: 1.90, 70: 1.87, 71: 1.87,
  72: 1.75, 73: 1.70, 74: 1.62, 75: 1.51, 76: 1.44, 77: 1.41, 78: 1.36,
  80: 1.32, 81: 1.45, 82: 1.46, 83: 1.48, 84: 1.40,
  86: 1.50,
  88: 2.21, 89: 2.15, 90: 2.06, 91: 2.00, 92: 1.96, 93: 1.90, 94: 1.87
};

export function get_covalent_radius(Z) {
  Z = Z | 0;
  if (Z <= 0) return 0.77;

  // 1) minimal fallback table
  const f = fallback_rcovalent(Z);
  if (f && f > 0) return f;

  // 2) small extra table
  const e = _RCOV_EXTRA[Z];
  if (e && e > 0) return e;

  // 3) heuristic by Z range (coarse but stable)
  if (Z <= 2) return 0.28;
  if (Z <= 10) return 0.70;
  if (Z <= 18) return 1.05;
  if (Z <= 36) return 1.25;
  if (Z <= 54) return 1.45;
  if (Z <= 86) return 1.65;
  return 1.80;
}

function _atomZ(a) {
  if (!a || typeof a !== 'object') return 0;
  if ('Z' in a) {
    const z = a.Z | 0;
    if (z > 0) return z;
  }
  if ('element' in a) return _symbol_to_Z(a.element);
  return 0;
}

// Spatial hash key (3 integers -> string). Works fine for view-sized subsets.
function _cellKey(ix, iy, iz) { return ix + ',' + iy + ',' + iz; }

/**
 * Guess bonds by distance with covalent radii heuristic.
 *
 * @param {Array<Object>} atoms - [{Z|element,x,y,z}, ...] in Å.
 * @param {Object} opts - { periodic?:boolean, cell?:any, maxNeighborsPerAtom?:number, tolA?:number, minDistA?:number }
 * @returns {Array<Array<number>>} bonds - [[i,j,1], ...] (single bonds only)
 */
export function guess_bonds_by_distance(atoms, opts) {
  opts = opts || {};
  if (!Array.isArray(atoms) || atoms.length < 2) return [];

  const tolA = (Number.isFinite(opts.tolA) ? opts.tolA : 0.45);
  const minDistA = (Number.isFinite(opts.minDistA) ? opts.minDistA : 0.63);
  const maxNeighbors = (Number.isFinite(opts.maxNeighborsPerAtom) ? (opts.maxNeighborsPerAtom | 0) : 12);

  const n = atoms.length;

  // Prepare radii + positions and choose a grid cell size.
  const R = new Float32Array(n);
  const X = new Float32Array(n);
  const Y = new Float32Array(n);
  const Zc = new Float32Array(n);

  let maxR = 0;
  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    const z = _atomZ(a);
    const r = get_covalent_radius(z);
    R[i] = r;
    if (r > maxR) maxR = r;

    X[i] = Number(a.x) || 0;
    Y[i] = Number(a.y) || 0;
    Zc[i] = Number(a.z) || 0;
  }

  // Grid cell size: large enough so that any potential neighbor is within +/-1 cells.
  const cellSize = Math.max(1.0, (2.0 * maxR + tolA));
  const inv = 1.0 / cellSize;

  // Build grid: Map<"ix,iy,iz", int[]>
  const grid = new Map();
  for (let i = 0; i < n; i++) {
    const ix = Math.floor(X[i] * inv);
    const iy = Math.floor(Y[i] * inv);
    const iz = Math.floor(Zc[i] * inv);
    const k = _cellKey(ix, iy, iz);
    let arr = grid.get(k);
    if (!arr) { arr = []; grid.set(k, arr); }
    arr.push(i);
  }

  const deg = (maxNeighbors > 0) ? new Uint16Array(n) : null;
  const bonds = [];

  // Early global cap to avoid runaway in dense crystals
  const maxBondsTotal = Math.max(1024, Math.min(2_000_000, (n * Math.max(1, maxNeighbors))));

  const minDist2 = minDistA * minDistA;

  for (let i = 0; i < n; i++) {
    if (deg && deg[i] >= maxNeighbors) continue;

    const ix0 = Math.floor(X[i] * inv);
    const iy0 = Math.floor(Y[i] * inv);
    const iz0 = Math.floor(Zc[i] * inv);

    for (let dz = -1; dz <= 1; dz++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dx = -1; dx <= 1; dx++) {
          const key = _cellKey(ix0 + dx, iy0 + dy, iz0 + dz);
          const list = grid.get(key);
          if (!list) continue;

          for (let t = 0; t < list.length; t++) {
            const j = list[t];
            if (j <= i) continue;
            if (deg && (deg[i] >= maxNeighbors || deg[j] >= maxNeighbors)) continue;

            const dxA = X[j] - X[i];
            const dyA = Y[j] - Y[i];
            const dzA = Zc[j] - Zc[i];
            const d2 = dxA * dxA + dyA * dyA + dzA * dzA;

            if (d2 < minDist2) continue;

            const cut = (R[i] + R[j] + tolA);
            const cut2 = cut * cut;
            if (d2 <= cut2) {
              bonds.push([i, j, 1]);
              if (deg) { deg[i]++; deg[j]++; }
              if (bonds.length >= maxBondsTotal) return bonds;
            }
          }
        }
      }
    }
  }

  return bonds;
}
