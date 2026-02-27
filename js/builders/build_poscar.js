// js/builders/build_poscar.js
// POSCAR/CONTCAR system builder used by phases.js.

import { load_poscar_system } from '../io/poscar_io.js';

export async function build_from_poscar(spec) {
  const sys = await load_poscar_system(spec.path);

  // POSCAR has no bonds => bonds=null (guess allowed later in provider.getView).
  const out = [sys.atoms || sys.atoms_cart || [], null, sys.title || 'POSCAR'];

  // Attach periodic extras for phases.js -> PeriodicProvider creation.
  out._tem_kind = 'periodic';
  if (Array.isArray(sys.cell)) out._tem_cell = sys.cell;
  if (Array.isArray(sys.atoms_cart)) out._tem_atoms_cart = sys.atoms_cart;
  if (Array.isArray(sys.atoms_frac)) out._tem_atoms_frac = sys.atoms_frac;
  if (sys.cell_from) out._tem_cell_from = sys.cell_from;
  return out;
}
