// js/builders/build_csv.js
// Strict CSV system builder used by phases.js.

import { load_csv_strict } from '../io/csv_io.js';

export async function build_from_csv(spec) {
  const pair = await load_csv_strict(spec.path);
  const atoms = pair[0];
  const bonds = pair[1];
  return [atoms, bonds, 'CSV'];
}
