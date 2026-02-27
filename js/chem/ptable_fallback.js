// js/chem/ptable_fallback.js
// Small, dependency-free fallbacks for CIF/CSV paths so they do not depend on RDKit WASM being ready.
// Kept intentionally minimal; extend only when needed.

export const FALLBACK_Z = {
  H: 1, He: 2,
  Li: 3, Be: 4, B: 5, C: 6, N: 7, O: 8, F: 9, Ne: 10,
  Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, Ar: 18,
  K: 19, Ca: 20,
  Fe: 26, Co: 27, Ni: 28, Cu: 29, Zn: 30,
  Br: 35, Ag: 47, I: 53, Au: 79
};

export const FALLBACK_RCOV = {
  1: 0.31, 2: 0.28,
  3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 10: 0.58,
  11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 18: 1.06,
  19: 2.03, 20: 1.76,
  26: 1.24, 27: 1.18, 28: 1.17, 29: 1.22, 30: 1.20,
  35: 1.20, 47: 1.45, 53: 1.39, 79: 1.36
};

export function norm_sym(sym) {
  sym = String(sym || '').trim();
  if (!sym) return '';
  const m = sym.match(/^([A-Za-z]{1,2})/);
  sym = m ? m[1] : sym;
  return sym.length === 1
    ? sym[0].toUpperCase()
    : sym[0].toUpperCase() + sym.slice(1).toLowerCase();
}

export function fallback_symbol_to_Z(sym) {
  const n = norm_sym(sym);
  return FALLBACK_Z[n] || 0;
}

export function fallback_rcovalent(Z) {
  return FALLBACK_RCOV[Z] || 0;
}
