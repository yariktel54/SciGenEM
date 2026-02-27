// js/builders/build_pdb.js
// PDB builder: converts PDB text into the internal System format (atoms + bonds semantics).

import { load_pdb } from '../io/pdb_io.js';
import { is_dev_mode } from '../system/system_contract.js';

// Bonds policy (strict):
//  - If CONECT exists => explicit bonds (array) => guessing is forbidden.
//  - If CONECT absent:
//      * large structures => bonds=[] (explicit none) to avoid accidental guessing
//      * small structures => bonds=null (unknown) to allow late fallback guessing (provider.getView)
const PDB_MAX_ATOMS_FOR_GUESS = 800; // "few hundred" safe range; tweak if needed

function title_from_path(path) {
  const s = String(path || '');
  if (!s) return 'PDB';
  // prefer fragment (blob:...#file.pdb)
  const frag = s.split('#').pop();
  const base = (frag || s).split('/').pop();
  return (base && base.trim()) ? base.trim() : 'PDB';
}

export async function build_from_pdb(spec) {
  spec = spec || {};
  const parsed = await load_pdb(spec.path);

  const atoms = Array.isArray(parsed.atoms) ? parsed.atoms : [];

  let bonds = null;
  if (parsed.hasConect) {
    bonds = Array.isArray(parsed.bonds) ? parsed.bonds : [];
  } else {
    bonds = (atoms.length <= PDB_MAX_ATOMS_FOR_GUESS) ? null : [];
  }

  const title = (parsed.title && String(parsed.title).trim()) ? String(parsed.title).trim() : title_from_path(spec.path);

  // DEV: one short line to verify altLoc/occupancy de-dup worked.
  if (is_dev_mode()) {
    const b = (bonds === null) ? 'null' : String((Array.isArray(bonds) ? bonds.length : 0));
    const before = parsed?.meta?.atomsBefore;
    const after = parsed?.meta?.atomsAfter;
    console.debug(`[TEM] PDB altLoc filter: atomsBefore=${before}; atomsAfter=${after}; bonds=${b}`);
  }

  // Attach lightweight warnings (kept optional; phases.js will merge these if present).
  const out = [atoms, bonds, title];
  out.meta = {
    warnings: [
      parsed.hasConect ? 'PDB: CONECT present (explicit bonds)' : (atoms.length > PDB_MAX_ATOMS_FOR_GUESS ? 'PDB: no CONECT (bonds=[])' : 'PDB: no CONECT (bonds=null => guess allowed)')
    ]
  };
  return out;
}
