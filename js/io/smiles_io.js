// js/io/smiles_io.js
// create_molecule_from_smiles(smiles), get_atoms_with_coords(mol)
import { RDKit } from "../chem/rdkit_wrap.js";
import { is_dev_mode } from "../system/system_contract.js";

function separate_fragments_inplace(atoms, bonds, gap = 8.0) {
  if (!atoms.length) return false;

  const g = new Map();
  for (let i = 0; i < atoms.length; i++) g.set(i, []);
  for (const [i, j] of bonds) {
    if (i >= 0 && j >= 0 && i < atoms.length && j < atoms.length) {
      g.get(i).push(j);
      g.get(j).push(i);
    }
  }

  const comps = [];
  const seen = new Array(atoms.length).fill(false);
  const stack = [];
  for (let s = 0; s < atoms.length; s++) {
    if (seen[s]) continue;
    stack.length = 0;
    stack.push(s);
    seen[s] = true;
    const comp = [];
    while (stack.length) {
      const v = stack.pop();
      comp.push(v);
      for (const u of g.get(v)) {
        if (!seen[u]) {
          seen[u] = true;
          stack.push(u);
        }
      }
    }
    comps.push(comp);
  }
  if (comps.length <= 1) return false;

  let shiftX = 0;
  for (const comp of comps) {
    const cx = comp.reduce((s, i) => s + atoms[i].x, 0) / comp.length;
    const cy = comp.reduce((s, i) => s + atoms[i].y, 0) / comp.length;
    for (const i of comp) {
      atoms[i].x = atoms[i].x - cx + shiftX;
      atoms[i].y = atoms[i].y - cy;
    }
    shiftX += gap;
  }
  return true;
}

function _tem_push_warning(mol, msg) {
  try {
    if (!mol) return;
    if (!Array.isArray(mol.__tem_warnings)) mol.__tem_warnings = [];
    const text = String(msg || "").trim();
    if (!text) return;
    if (mol.__tem_warnings.indexOf(text) >= 0) return;
    mol.__tem_warnings.push(text);
  } catch (_) {}
}

function _tem_dev_enabled() {
  try {
    return !!(
      typeof window !== "undefined" &&
      window &&
      window.TEM_DEV === true
    );
  } catch (_) {
    return false;
  }
}

function _smiles_io_dev() {
  return is_dev_mode() || _tem_dev_enabled();
}

function _smiles_io_log(...args) {
  if (_smiles_io_dev()) console.log("[SMILES IO]", ...args);
}

function _safe_center_xy(atoms) {
  if (!Array.isArray(atoms) || !atoms.length) return;
  const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length;
  const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length;
  for (const a of atoms) {
    a.x -= cx;
    a.y -= cy;
  }
}

function _clone_atoms(atoms) {
  return Array.isArray(atoms) ? atoms.map((a) => ({ ...a })) : [];
}

function _clone_bonds(bonds) {
  return Array.isArray(bonds)
    ? bonds.map((b) => (Array.isArray(b) ? b.slice(0, 3) : b))
    : [];
}

function _normalize_smiles_io_options(use_2d_if_fail, options) {
  let cfg = { use_2d_if_fail: true, strictChemistry: false };
  if (typeof use_2d_if_fail === "object" && use_2d_if_fail !== null) {
    options = use_2d_if_fail;
    use_2d_if_fail = undefined;
  }
  if (typeof use_2d_if_fail === "boolean") cfg.use_2d_if_fail = use_2d_if_fail;
  if (options && typeof options === "object") {
    if (typeof options.use_2d_if_fail === "boolean")
      cfg.use_2d_if_fail = options.use_2d_if_fail;
    if (options.strictChemistry === true) cfg.strictChemistry = true;
  }
  return cfg;
}

function _make_strict_error(trust) {
  const reasons =
    trust && Array.isArray(trust.rejectReasons) ? trust.rejectReasons : [];
  const msg =
    "SMILES strict chemistry gate rejected: " +
    (reasons.length ? reasons.join(", ") : "strict_rejected");
  const err = new Error(msg);
  err.rejectReasons = reasons.slice(0);
  err.trust = trust || null;
  return err;
}

function _attach_smiles_backend_to_mol(mol, backend, geometry, trust) {
  if (!mol) return;
  try {
    mol.__tem_smiles_backend = backend || null;
    mol.__tem_smiles_geometry = geometry || null;
    mol.__tem_smiles_trust = trust || null;
    if (!Array.isArray(mol.__tem_warnings)) mol.__tem_warnings = [];
    const warnings = [];
    if (backend && Array.isArray(backend.warnings))
      warnings.push(...backend.warnings);
    if (geometry && Array.isArray(geometry.warnings))
      warnings.push(...geometry.warnings);
    if (trust && Array.isArray(trust.rejectReasons))
      warnings.push(...trust.rejectReasons);
    for (const msg of warnings) _tem_push_warning(mol, msg);
  } catch (_) {}
}

export async function create_molecule_from_smiles(smiles, options) {
  const opts = _normalize_smiles_io_options(undefined, options);
  const DEV = _smiles_io_dev();
  const dwarn = (...a) => {
    if (DEV) console.warn(...a);
  };

  _smiles_io_log("input =", smiles);

  if (window.RDKitReady && typeof window.RDKitReady.then === "function") {
    if (!window.RDKit) window.RDKit = RDKit;
    if (typeof RDKit.ensure_smiles_backends_ready === "function") {
      const ready = await RDKit.ensure_smiles_backends_ready();
      if (ready && ready.rdkit && typeof RDKit.setModule === "function") {
        RDKit.setModule(ready.rdkit);
      }
    } else {
      const mod = await window.RDKitReady;
      if (typeof RDKit.setModule === "function") RDKit.setModule(mod);
    }
  } else {
    dwarn("RDKitReady не знайдено або ще не завантажений.");
    throw new Error("RDKit не готовий");
  }

  const backend = RDKit.build_smiles_backend_result(smiles);
  const geometry = RDKit.finalize_smiles_geometry(backend);
  const trust = RDKit.assess_smiles_trust({ backend, geometry });
  const mol = backend && backend.mol ? backend.mol : null;

  const replacementStatus =
    typeof RDKit.replacement_get_status === "function"
      ? RDKit.replacement_get_status()
      : null;

  _smiles_io_log("replacement status", replacementStatus);

  _smiles_io_log("backend summary", {
    ok: !!backend?.ok,
    atomCountFromMol: backend?.diagnostics?.atomCountFromMol ?? 0,
    atomCountFromAtoms: backend?.diagnostics?.atomCountFromAtoms ?? 0,
    bondCountFromMol: backend?.diagnostics?.bondCountFromMol ?? 0,
    bondCountFromBonds: backend?.diagnostics?.bondCountFromBonds ?? 0,
    coordsAvailable: !!backend?.stages?.coordsAvailable,
    hydrogenMode: backend?.identity?.hydrogenMode || "unsupported",
    explicitHydrogenCount: backend?.identity?.explicitHydrogenCount ?? 0,
    implicitHydrogenCount: backend?.identity?.implicitHydrogenCount ?? 0,
    uniqueSymbols: backend?.identity?.uniqueSymbols || [],
    uiUniqueSymbols: backend?.identity?.uiUniqueSymbols || [],
    authoritativeTopologySource:
      backend?.authoritativeTopologySource || "fallback",
    authoritativeHydrogenSource:
      backend?.authoritativeHydrogenSource || "fallback",
    authoritativeGeometrySource:
      backend?.authoritativeGeometrySource || "fallback",
    authoritativeMolKind: backend?.authoritativeMolKind || "fallback",
    warnings: backend?.warnings || [],
    diagnostics: backend?.diagnostics || null,
  });

  _smiles_io_log("geometry summary", {
    ok: !!geometry?.ok,
    geometryMode: geometry?.geometryMode || "failed",
    coordsAvailable: !!geometry?.coordsAvailable,
    degenerateCoords: !!geometry?.degenerateCoords,
    usedFallbackCoords: !!geometry?.usedFallbackCoords,
    warnings: geometry?.warnings || [],
  });

  _smiles_io_log("trust summary", {
    strictOk: !!trust?.strictOk,
    trustLevel: trust?.trustLevel || "rejected",
    rejectReasons: trust?.rejectReasons || [],
  });

  if (!backend || !backend.ok || !mol) {
    throw new Error("Не вдалося зчитати SMILES: " + smiles);
  }

  _attach_smiles_backend_to_mol(mol, backend, geometry, trust);

  if (opts.strictChemistry && (!trust || !trust.strictOk)) {
    throw _make_strict_error(trust);
  }

  return mol;
}

export function get_atoms_with_coords(mol, use_2d_if_fail = true, options) {
  const opts = _normalize_smiles_io_options(use_2d_if_fail, options);

  let backend = null;
  let geometry = null;
  let trust = null;

  try {
    backend = mol && mol.__tem_smiles_backend ? mol.__tem_smiles_backend : null;
    geometry =
      mol && mol.__tem_smiles_geometry ? mol.__tem_smiles_geometry : null;
    trust = mol && mol.__tem_smiles_trust ? mol.__tem_smiles_trust : null;
  } catch (_) {
    backend = null;
    geometry = null;
    trust = null;
  }

  if (!geometry || !geometry.ok) {
    geometry = RDKit.finalize_smiles_geometry(backend || mol);
    trust = RDKit.assess_smiles_trust({ backend, geometry });
    _attach_smiles_backend_to_mol(mol, backend, geometry, trust);
  }

  if (opts.strictChemistry && (!trust || !trust.strictOk)) {
    throw _make_strict_error(trust);
  }

  let atoms = _clone_atoms(geometry && geometry.atoms);
  let bonds = _clone_bonds(geometry && geometry.bonds);
  const usedFallbackCoords = !!(geometry && geometry.usedFallbackCoords);

  if (!atoms.length) {
    throw new Error("❌ Молекула не має координат.");
  }

  if (Array.isArray(bonds) && bonds.length < atoms.length - 1) {
    if (separate_fragments_inplace(atoms, bonds, 8.0)) {
      _tem_push_warning(mol, "fragment_separation_applied");
    }
  }

  _safe_center_xy(atoms);

  const uniq = Array.from(
    new Set(
      atoms
        .map((a) => (Number.isFinite(a?.Z) ? a.Z | 0 : 0))
        .filter((z) => z > 0),
    ),
  );

  _smiles_io_log("final atoms/bonds summary", {
    atomCount: atoms.length,
    bondCount: Array.isArray(bonds) ? bonds.length : 0,
    uniqueZ: uniq,
    hydrogenCount: atoms.filter((a) => (a?.Z | 0) === 1).length,
    usedFallbackCoords,
    geometryMode: geometry?.geometryMode || "failed",
    strictOk: !!(trust && trust.strictOk),
  });

  return [atoms, Array.isArray(bonds) ? bonds : []];
}
