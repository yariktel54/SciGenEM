// js/builders/build_smiles.js
// SMILES system builder (RDKit only) used by phases.js.

import {
  create_molecule_from_smiles,
  get_atoms_with_coords,
} from "../io/smiles_io.js";

const STRICT_SMILES_CHEMISTRY = true;

export async function build_from_smiles(spec) {
  const rawSmiles = String((spec && spec.smiles) || "");
  const inputSmiles = rawSmiles.trim();
  const buildSmiles = inputSmiles || "O";

  let mol;
  try {
    mol = await create_molecule_from_smiles(buildSmiles, {
      strictChemistry: STRICT_SMILES_CHEMISTRY,
    });
  } catch (e) {
    const err = new Error(
      `SMILES build failed: ${inputSmiles || "<empty input>"}`,
    );
    err.cause = e;
    err.rejectReasons = Array.isArray(e && e.rejectReasons)
      ? e.rejectReasons.slice(0)
      : [];
    err.trust = e && e.trust ? e.trust : null;
    throw err;
  }

  const pair = get_atoms_with_coords(mol, {
    strictChemistry: STRICT_SMILES_CHEMISTRY,
  });
  const atoms = pair[0];
  const bonds = pair[1];

  let titleSmiles = buildSmiles;
  try {
    const used = mol?.__tem_smiles_backend?.smiles;
    if (typeof used === "string" && used.trim()) titleSmiles = used.trim();
  } catch (_) {}

  const out = [atoms, bonds, `SMILES: ${titleSmiles || "O"}`];
  try {
    const w = mol?.__tem_warnings;
    const trust = mol?.__tem_smiles_trust;
    const meta = {};
    if (Array.isArray(w) && w.length) meta.warnings = w.slice(0, 16);
    if (trust && typeof trust === "object") {
      meta.trust = {
        strictOk: !!trust.strictOk,
        trustLevel: trust.trustLevel || "rejected",
        rejectReasons: Array.isArray(trust.rejectReasons)
          ? trust.rejectReasons.slice(0, 16)
          : [],
      };
    }
    if (Object.keys(meta).length) out.meta = meta;
  } catch (_) {}
  return out;
}
