// js/builders/build_xyz.js
// XYZ system builder used by phases.js.

import { load_xyz } from '../io/xyz_io.js';

export async function build_from_xyz(spec) {
  const sys = await load_xyz(spec.path);
  // XYZ has no bonds => bonds=null (guess allowed later).
  return [sys.atoms, sys.bonds, sys.title || 'XYZ'];
}
