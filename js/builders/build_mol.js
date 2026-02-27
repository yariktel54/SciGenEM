// js/builders/build_mol.js
// MOL (Molfile) system builder used by phases.js.

import { load_mol } from '../io/mol_io.js';

export async function build_from_mol(spec) {
  const sys = await load_mol(spec.path);
  // MOL has explicit bonds -> bonds is an array (guess forbidden later).
  return [sys.atoms, sys.bonds, sys.title || 'MOL'];
}
