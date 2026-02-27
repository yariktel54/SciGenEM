// js/builders/build_smiles.js
// SMILES system builder (RDKit only) used by phases.js.

import { create_molecule_from_smiles, get_atoms_with_coords } from '../io/smiles_io.js';

export async function build_from_smiles(spec) {
  let sm = (spec.smiles || '').trim();
  if (!sm) sm = 'O';

  let mol;
  try {
    mol = await create_molecule_from_smiles(sm);
  } catch (e) {
    sm = 'O';
    mol = await create_molecule_from_smiles(sm);
    try {
      if (!Array.isArray(mol.__tem_warnings)) mol.__tem_warnings = [];
      mol.__tem_warnings.push('SMILES input failed; fallback to O');
    } catch (_) {}
  }

  // Bonds come strictly from RDKit here.
  const pair = get_atoms_with_coords(mol);
  const atoms = pair[0];
  const bonds = pair[1];

    const out = [atoms, bonds, `SMILES: ${sm}`];
  // Pass through any non-fatal warnings collected during SMILES parsing / coord generation.
  try {
    const w = mol?.__tem_warnings;
    if (Array.isArray(w) && w.length) out.meta = { warnings: w.slice(0, 8) };
  } catch (_) {}
  return out;
}
