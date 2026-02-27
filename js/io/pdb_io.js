// js/io/pdb_io.js
// Minimal, robust PDB (v3.x-style) loader.
// Scope:
//  - Parse ATOM/HETATM coordinates (fixed-width columns)
//  - Handle altLoc/occupancy de-dup (avoid double-darkening)
//  - Support MODEL/ENDMDL: take only MODEL 1 (or first model)
//  - Parse CONECT into explicit single bonds (type=1)
//
// IMPORTANT bonds semantics are decided in the builder (build_pdb.js):
//  - If CONECT exists -> bonds is an array (explicit topology)
//  - If CONECT absent -> builder sets bonds = [] for large systems, or null for small (guess allowed)

import { norm_sym, fallback_symbol_to_Z } from '../chem/ptable_fallback.js';

// Full periodic table mapping (symbol -> Z) as a robust fallback (no RDKit needed).
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

function symbol_to_Z(sym) {
  const n = norm_sym(sym);
  return fallback_symbol_to_Z(n) || _Z_BY_SYM[n] || 0;
}

function _slice(line, a, b) {
  if (!line) return '';
  if (a >= line.length) return '';
  return line.slice(a, Math.min(b, line.length));
}

function _toInt(s) {
  const v = parseInt(String(s || '').trim(), 10);
  return Number.isFinite(v) ? (v | 0) : 0;
}

function _toFloat(s, def = 0.0) {
  const v = parseFloat(String(s || '').trim());
  return Number.isFinite(v) ? v : def;
}

function _alt_priority(altLoc) {
  // Policy:
  //  a) keep altLoc=' ' OR 'A' (priority)
  //  b) else pick max occupancy
  // Encode as priority score: ' ' highest, then 'A', then others.
  if (altLoc === ' ' || altLoc === '' || altLoc == null) return 2;
  if (altLoc === 'A') return 1;
  return 0;
}

function _infer_element_from_atom_name(atomName4) {
  // atomName4: columns 13-16, includes leading spaces.
  const raw = String(atomName4 || '');
  if (!raw) return '';

  // PDB convention:
  //  - one-letter elements are right-justified => raw[0] is space, element at raw[1]
  //  - two-letter elements are left-justified => raw[0] is a letter
  // This avoids the common CA (C-alpha) vs Ca (Calcium) pitfall.
  if (raw.length >= 2 && raw[0] === ' ' && /[A-Za-z]/.test(raw[1])) {
    return norm_sym(raw[1]);
  }

  let t = raw.trim();
  if (!t) return '';

  // Strip leading digit (e.g., 1HG1 -> HG1)
  if (t.length >= 2 && /\d/.test(t[0])) t = t.slice(1);
  if (!t) return '';

  const c0 = t[0];
  const c1 = (t.length >= 2) ? t[1] : '';

  // If second char is a letter => candidate two-letter element, else one-letter.
  if (c1 && /[A-Za-z]/.test(c1)) return norm_sym(c0 + c1);
  return norm_sym(c0);
}

function _parse_conect_fields(line) {
  // CONECT record uses 5-char integer fields after column 6.
  // Columns (1-based):
  //  7-11: atom serial
  // 12-16,17-21,22-26,27-31: bonded atom serials (may extend in practice)
  const out = [];
  const src = _toInt(_slice(line, 6, 11));
  if (!src) return out;
  for (let pos = 11; pos + 5 <= line.length; pos += 5) {
    const v = _toInt(_slice(line, pos, pos + 5));
    if (v) out.push([src, v]);
  }
  return out;
}

export function parse_pdb_text(txt, opts) {
  opts = opts || {};
  const text = String(txt || '').replace(/^\uFEFF/, '');
  const lines = text.split(/\r?\n/);

  // MODEL handling
  let hasModel = false;
  let selectedModel = null; // numeric or null
  let currentModel = null;
  let inModel = false;
  let afterSelectedEnd_beforeNextModel = false;

  // Title
  const titleParts = [];
  let headerLine = '';

  // Atoms with altLoc de-dup
  const keyToIndex = new Map();
  const atomsSel = [];
  let atomLinesTotal = 0;

  // CONECT
  const conectPairs = [];
  let hasConect = false;

  function acceptAtomsNow() {
    if (!hasModel) return true;
    return inModel && (currentModel === selectedModel);
  }

  function acceptConectNow() {
    if (!hasModel) return true;
    if (inModel && (currentModel === selectedModel)) return true;
    if (!inModel && afterSelectedEnd_beforeNextModel) return true;
    return false;
  }

  for (let li = 0; li < lines.length; li++) {
    const line = lines[li];
    if (!line) continue;

    // Record name is columns 1-6
    const rec6 = line.length >= 6 ? line.slice(0, 6) : line;

    if (rec6 === 'MODEL ') {
      hasModel = true;
      // Model serial is columns 11-14 (10..14)
      let mn = _toInt(_slice(line, 10, 14));
      if (!mn) mn = (selectedModel == null) ? 1 : mn;

      if (selectedModel == null) {
        selectedModel = mn;
        afterSelectedEnd_beforeNextModel = false;
      } else if (afterSelectedEnd_beforeNextModel && mn !== selectedModel) {
        // Any later model starts => stop accepting CONECT outside the selected model
        afterSelectedEnd_beforeNextModel = false;
      }

      currentModel = mn;
      inModel = true;
      continue;
    }

    if (rec6 === 'ENDMDL') {
      if (hasModel && inModel && currentModel === selectedModel) {
        // After first (selected) model ends, some files put CONECT right after ENDMDL.
        afterSelectedEnd_beforeNextModel = true;
      }
      inModel = false;
      currentModel = null;
      continue;
    }

    if (rec6 === 'TITLE ') {
      // columns 11-80 are title text
      const t = _slice(line, 10, 80).trim();
      if (t) titleParts.push(t);
      continue;
    }

    if (rec6 === 'HEADER') {
      // Keep as fallback if TITLE is absent.
      const h = _slice(line, 10, 80).trim();
      if (h) headerLine = h;
      continue;
    }

    if (rec6 === 'ATOM  ' || rec6 === 'HETATM') {
      if (!acceptAtomsNow()) continue;
      atomLinesTotal++;

      // Fixed-width fields (PDB v3.x)
      const serial = _toInt(_slice(line, 6, 11));
      const atomName4 = _slice(line, 12, 16); // keep spacing for element inference
      const atomName = atomName4.trim();
      const altLoc = (line.length > 16) ? line[16] : ' ';
      const resName = _slice(line, 17, 20).trim();
      const chainID = (line.length > 21) ? String(line[21]) : '';
      const resSeqStr = _slice(line, 22, 26).trim();
      const resSeq = _toInt(resSeqStr);
      const iCode = (line.length > 26) ? String(line[26]).trim() : '';

      const x = _toFloat(_slice(line, 30, 38), 0.0);
      const y = _toFloat(_slice(line, 38, 46), 0.0);
      const z = _toFloat(_slice(line, 46, 54), 0.0);
      const occupancy = _toFloat(_slice(line, 54, 60), 0.0);
      const tempFactor = _toFloat(_slice(line, 60, 66), 0.0);

      let el = _slice(line, 76, 78).trim();
      el = norm_sym(el);
      if (!el) el = _infer_element_from_atom_name(atomName4);

      const Z = symbol_to_Z(el);

      // Grouping key for altLoc de-dup
      const key = chainID + '|' + resSeqStr + '|' + iCode + '|' + resName + '|' + atomName;
      const p = _alt_priority(altLoc);

      const cand = {
        serial,
        Z,
        element: el,
        x, y, z,
        meta: {
          chain: chainID,
          resName,
          resSeq,
          iCode,
          atomName,
          occupancy,
          tempFactor,
          altLoc: (altLoc === ' ' ? '' : altLoc)
        },
        _occ: occupancy,
        _altP: p
      };

      const prevIdx = keyToIndex.get(key);
      if (prevIdx == null) {
        keyToIndex.set(key, atomsSel.length);
        atomsSel.push(cand);
      } else {
        const prev = atomsSel[prevIdx];
        const better = (cand._altP > prev._altP) || (cand._altP === prev._altP && cand._occ > prev._occ);
        if (better) atomsSel[prevIdx] = cand;
      }

      continue;
    }

    if (rec6 === 'CONECT') {
      if (!acceptConectNow()) continue;
      hasConect = true;
      const pairs = _parse_conect_fields(line);
      if (pairs && pairs.length) conectPairs.push(...pairs);
      continue;
    }

    // If a later MODEL begins after we finished MODEL 1, stop accepting CONECT outside it.
    if (hasModel && afterSelectedEnd_beforeNextModel && rec6 === 'MODEL ') {
      afterSelectedEnd_beforeNextModel = false;
    }
  }

  // Strip helper fields
  const atoms = new Array(atomsSel.length);
  for (let i = 0; i < atomsSel.length; i++) {
    const a = atomsSel[i];
    atoms[i] = { Z: a.Z | 0, x: +a.x || 0, y: +a.y || 0, z: +a.z || 0, element: a.element, meta: a.meta };
  }

  // Build bonds from CONECT (explicit) with serial->index mapping.
  let bonds = [];
  if (hasConect && conectPairs.length) {
    const serialToIdx = new Map();
    for (let i = 0; i < atomsSel.length; i++) {
      const s = atomsSel[i].serial | 0;
      if (s) serialToIdx.set(s, i);
    }

    const seen = new Set();
    const out = [];

    for (let k = 0; k < conectPairs.length; k++) {
      const s1 = conectPairs[k][0] | 0;
      const s2 = conectPairs[k][1] | 0;
      const i = serialToIdx.get(s1);
      const j = serialToIdx.get(s2);
      if (i == null || j == null) continue;
      if (i === j) continue;
      const a = i < j ? i : j;
      const b = i < j ? j : i;
      const key = a + ',' + b;
      if (seen.has(key)) continue;
      seen.add(key);
      out.push([a, b, 1]);
    }

    bonds = out;
  }

  const title = (titleParts.length ? titleParts.join(' ') : headerLine).trim();

  return {
    title,
    atoms,
    hasConect: !!hasConect,
    bonds,
    meta: {
      atomsBefore: atomLinesTotal,
      atomsAfter: atoms.length,
      model: (selectedModel == null ? null : selectedModel)
    }
  };
}

export async function load_pdb(pathOrText) {
  const s = String(pathOrText || '').trim();
  let txt = s;

  if (/^(https?:\/\/|blob:)/i.test(s)) {
    const r = await fetch(s.split('#')[0]);
    if (!r || !r.ok) throw new Error('PDB: fetch failed');
    txt = await r.text();
  }

  return parse_pdb_text(txt);
}
