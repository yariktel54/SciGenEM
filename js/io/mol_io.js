// js/io/mol_io.js
// Molfile (.mol) loader via RDKit WASM (no full V2000/V3000 parser).
//
// Rules:
//  - MOL has explicit bonds -> bonds MUST be an ARRAY (guess forbidden).
//  - We rely on RDKit to read molblock and (if needed) convert to a canonical molblock.
//  - Extraction is done from RDKit-exported molblock (V2000 layout), so we don't implement V3000 parsing.

import { RDKit } from '../chem/rdkit_wrap.js';
import { is_dev_mode } from '../system/system_contract.js';
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

function looksLikeRawMolText(s) {
  s = String(s || '');
  if (!s.includes('\n')) return false;
  if (/(^|\n)M\s{2}END(\s|$)/.test(s)) return true;
  const lines = s.split(/\r?\n/);
  if (lines.length >= 4) {
    const counts = String(lines[3] || '');
    if (/V2000|V3000/i.test(counts)) return true;
  }
  return false;
}

function guessTitleFromPathOrText(pathOrText, rawText) {
  try {
    const lines = String(rawText || '').split(/\r?\n/);
    const t0 = String(lines[0] || '').trim();
    if (t0) return t0;
  } catch (_) {}
  const s = String(pathOrText || '');
  if (!s) return 'MOL';
  const frag = s.split('#').pop() || s;
  const base = frag.split('/').pop() || frag;
  return base || 'MOL';
}

async function ensureRdkitReady() {
  if (window.RDKitReady && typeof window.RDKitReady.then === 'function') {
    const mod = await window.RDKitReady;
    try { RDKit.setModule(mod); } catch (_) {}
    try {
      if (window.RDKit && typeof window.RDKit.setModule === 'function') window.RDKit.setModule(mod);
    } catch (_) {}
    return mod;
  }
  throw new Error('MOL: RDKit не готовий (RDKitReady відсутній)');
}

function createMolFromMolblock(mod, molblock) {
  let mol = null;

  // Many RDKit_minimal builds can autodetect molblock in get_mol().
  try { mol = mod.get_mol?.(molblock); } catch (_) { mol = null; }

  // Try with options (some builds accept JSON opts)
  if (!mol) {
    try { mol = mod.get_mol?.(molblock, JSON.stringify({ sanitize: true })); } catch (_) { mol = null; }
  }

  // Try constructor (some builds accept molblock)
  if (!mol && mod.Mol) {
    try { mol = new mod.Mol(molblock); } catch (_) { mol = null; }
  }

  try {
    if (mol && typeof mol.is_valid === 'function' && !mol.is_valid()) mol = null;
  } catch (_) {}

  return mol;
}

function parse_atoms_bonds_from_v2000_molblock(mb) {
  const lines = String(mb || '').split(/\r?\n/);
  if (lines.length < 4) return { atoms: [], bonds: [] };

  // Counts line (row 4, index 3)
  const counts = lines[3] || '';
  const nat = parseInt(counts.slice(0, 3), 10) || 0;
  const nb  = parseInt(counts.slice(3, 6), 10) || 0;
  if (nat <= 0 || lines.length < 4 + nat) return { atoms: [], bonds: [] };

  const atoms = [];
  for (let i = 0; i < nat; i++) {
    const ln = lines[4 + i] || '';
    const x = parseFloat(ln.slice(0, 10));
    const y = parseFloat(ln.slice(10, 20));
    const z = parseFloat(ln.slice(20, 30));
    const sym = ln.slice(31, 34).trim();
    const Z = symbol_to_Z(sym) || 6;
    atoms.push({ Z: Z | 0, x: Number.isFinite(x) ? +x : 0, y: Number.isFinite(y) ? +y : 0, z: Number.isFinite(z) ? +z : 0 });
  }

  const bonds = [];
  const bondStart = 4 + nat;
  for (let k = 0; k < nb; k++) {
    const ln = lines[bondStart + k] || '';
    const i = (parseInt(ln.slice(0, 3), 10) || 0) - 1;
    const j = (parseInt(ln.slice(3, 6), 10) || 0) - 1;
    const order = parseInt(ln.slice(6, 9), 10) || 1;
    if (i >= 0 && j >= 0 && i < nat && j < nat && i !== j) bonds.push([i, j, order]);
  }

  return { atoms, bonds };
}

export async function load_mol(pathOrText) {
  const DEV = is_dev_mode();
  const dlog = (...a) => { if (DEV) console.debug(...a); };

  let text = pathOrText ?? '';
  const s = String(text);

  // If it's a URL/blob, fetch it unless it already looks like raw mol text.
  if (!looksLikeRawMolText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
    const r = await fetch(s);
    if (!r.ok) throw new Error('MOL: не вдалося завантажити (' + r.status + ')');
    text = await r.text();
  }

  const raw = String(text || '');
  if (!raw.trim()) throw new Error('MOL: порожній вміст.');

  const mod = await ensureRdkitReady();
  const mol = createMolFromMolblock(mod, raw);
  if (!mol) throw new Error('MOL: RDKit не зміг зчитати molfile.');

  const title = guessTitleFromPathOrText(pathOrText, raw);

  // Canonical molblock from RDKit (typically V2000); we parse atoms+bonds from it.
  const mb = mol.get_molblock?.() ?? mol.MolToMolBlock?.() ?? '';
  const { atoms, bonds } = parse_atoms_bonds_from_v2000_molblock(mb);

  if (!atoms.length) {
    throw new Error('MOL: не вдалося отримати атоми з molblock (перевірте файл).');
  }

  // bonds are explicit in MOL => return array (even if empty)
  if (DEV) dlog(`[TEM] load mol: atoms=${atoms.length}; bonds=${bonds.length}; title=${title}`);

  return {
    title,
    kind: 'molecule',
    atoms,
    bonds
  };
}
