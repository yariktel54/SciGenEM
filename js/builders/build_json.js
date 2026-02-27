// js/builders/build_json.js
// JSON system builder used by phases.js.

import { load_json_system } from '../io/json_io.js';

export async function build_from_json(spec) {
  const sys = await load_json_system(spec.path);
  const out = [sys.atoms, sys.bonds, sys.title || 'JSON'];
  // CIF-like extras (optional): allow phases/provider to treat periodic JSON like CIF.
  if (Array.isArray(sys.atoms_cart)) out._tem_atoms_cart = sys.atoms_cart;
  if (Array.isArray(sys.atoms_frac)) out._tem_atoms_frac = sys.atoms_frac;
  if (sys.cell_from) out._tem_cell_from = sys.cell_from;

  // Attach small extras for phases.js (keeps legacy [atoms,bonds,title] contract intact).
  out._tem_kind = sys.kind || 'molecule';
  if (sys.cell) out._tem_cell = sys.cell;
  if (sys.cell_vectors) out._tem_cell_vectors = sys.cell_vectors;
  return out;
}
