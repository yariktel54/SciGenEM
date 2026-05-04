import { ensureOpenChemLibLoaded, getOpenChemLibModule, getOpenChemLibLoaderStatus } from "./openchemlib_loader.js";

// js/chem/rdkit_wrap.js
// Єдиний backend-adapter для RDKit_minimal.js:
// - raw layer: тонкі безпечні виклики до RDKit WASM
// - safe adapter layer: нормалізація, molblock parsing, diagnostics, fallbacks без хибної хімії

let _RDKit = null;

const ELEMENT_SYMBOLS = [
  null,
  "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
];

const SYMBOL_TO_Z = (() => {
  const out = Object.create(null);
  for (let z = 1; z < ELEMENT_SYMBOLS.length; z++) {
    const sym = ELEMENT_SYMBOLS[z];
    if (sym) out[sym] = z;
  }
  return out;
})();


const REPLACEMENT_BACKEND_NAME = "openchemlib";

let _replacementChem = {
  module: null,
  ready: null,
  source: "none",
  lastError: null
};

async function ensure_replacement_backend_loading() {
  if (_replacementChem.module) return _replacementChem.module;
  if (_replacementChem.ready) return _replacementChem.ready;

  _replacementChem.ready = Promise.resolve().then(async () => {
    const mod = await ensureOpenChemLibLoaded();
    _replacementChem.module = mod;
    const status = getOpenChemLibLoaderStatus();
    _replacementChem.source = status && status.source ? status.source : REPLACEMENT_BACKEND_NAME;
    _replacementChem.lastError = null;
    return mod;
  }).catch((e) => {
    _replacementChem.lastError = e;
    _replacementChem.ready = null;
    throw e;
  });

  return _replacementChem.ready;
}

function replacement_start_loading() {
  if (_replacementChem.module) return Promise.resolve(_replacementChem.module);
  if (_replacementChem.ready) return _replacementChem.ready;
  return ensure_replacement_backend_loading();
}

function replacement_get_module() {
  if (_replacementChem.module) return _replacementChem.module;
  const mod = getOpenChemLibModule();
  if (mod) {
    _replacementChem.module = mod;
    const status = getOpenChemLibLoaderStatus();
    if (status && status.source) _replacementChem.source = status.source;
    if (status && status.error) _replacementChem.lastError = new Error(status.error);
  }
  return _replacementChem.module;
}

function replacement_ready() {
  return !!replacement_get_module();
}

function replacement_get_status() {
  const loaderStatus = getOpenChemLibLoaderStatus();
  return {
    ready: !!(loaderStatus && loaderStatus.ready) || replacement_ready(),
    source:
      (loaderStatus && loaderStatus.source) ||
      _replacementChem.source ||
      REPLACEMENT_BACKEND_NAME,
    error:
      (loaderStatus && loaderStatus.error) ||
      (_replacementChem.lastError ? String(_replacementChem.lastError.message || _replacementChem.lastError) : "")
  };
}


async function ensure_smiles_backends_ready(options = {}) {
  let mod = _RDKit;
  if (!mod) {
    if (!(window.RDKitReady && typeof window.RDKitReady.then === "function")) {
      throw new Error("RDKitReady not available");
    }
    mod = await window.RDKitReady;
    _RDKit = mod;
  }

  const retryCount = Math.max(0, parse_int_safe(options.retryCount, 1));
  let lastError = null;
  for (let attempt = 0; attempt <= retryCount; attempt++) {
    try {
      await ensure_replacement_backend_loading();
      const status = replacement_get_status();
      if (status && status.ready) {
        return {
          rdkit: mod,
          replacementStatus: status
        };
      }
      lastError = new Error("SMILES replacement backend not ready");
    } catch (e) {
      lastError = e instanceof Error ? e : new Error(String(e));
    }
    if (attempt < retryCount) await Promise.resolve();
  }

  throw lastError || new Error("SMILES replacement backend not ready");
}

function replacement_is_wrapper_mol(mol) {
  return !!(mol && mol.__tem_adapter_wrapper === true);
}

function make_adapter_mol_wrapper(smiles) {
  return {
    __tem_adapter_wrapper: true,
    __tem_adapter_wrapper_type: "smiles_backend",
    __tem_input_smiles: String(smiles == null ? "" : smiles),
    delete() {}
  };
}


function replacement_clone_payload(payload) {
  const src = payload && typeof payload === "object" ? payload : null;
  if (!src) return null;
  return {
    ok: !!src.ok,
    backendName: src.backendName || REPLACEMENT_BACKEND_NAME,
    source: src.source || "replacement_backend",
    geometryOrigin: src.geometryOrigin || "replacement_2d",
    authoritativeTopologySource: src.authoritativeTopologySource || "fallback",
    authoritativeHydrogenSource: src.authoritativeHydrogenSource || "fallback",
    authoritativeGeometrySource: src.authoritativeGeometrySource || "fallback",
    authoritativeMolKind: src.authoritativeMolKind || "fallback",
    explicitHydrogenCount: Number.isFinite(src.explicitHydrogenCount) ? (src.explicitHydrogenCount | 0) : 0,
    implicitHydrogenCount: Number.isFinite(src.implicitHydrogenCount) ? (src.implicitHydrogenCount | 0) : 0,
    baseAtomCount: Number.isFinite(src.baseAtomCount) ? (src.baseAtomCount | 0) : 0,
    baseBondCount: Number.isFinite(src.baseBondCount) ? (src.baseBondCount | 0) : 0,
    hydrogenAddedCount: Number.isFinite(src.hydrogenAddedCount) ? (src.hydrogenAddedCount | 0) : 0,
    expandWorked: !!src.expandWorked,
    usedWrapperMol: !!src.usedWrapperMol,
    warnings: Array.isArray(src.warnings) ? src.warnings.slice(0) : [],
    diagnostics: src.diagnostics && typeof src.diagnostics === "object" ? { ...src.diagnostics } : {},
    identity: clone_identity_info(src.identity),
    atoms: Array.isArray(src.atoms) ? src.atoms.map(clone_atom_record) : [],
    bonds: Array.isArray(src.bonds) ? src.bonds.map(clone_bond_record).filter(Boolean) : [],
    coords: Array.isArray(src.coords) ? src.coords.map(clone_coord_record) : []
  };
}


function replacement_attach_to_mol(mol, payload) {
  if (!mol || !payload || !payload.ok) return;
  try {
    mol.__tem_replacement_topology = replacement_clone_payload(payload);
  } catch (_) {}
}

function replacement_get_cache(mol) {
  const cache = mol && mol.__tem_replacement_topology;
  if (cache && cache.ok) return cache;
  return null;
}

function replacement_count_h_atoms(records) {
  const arr = Array.isArray(records) ? records : [];
  let n = 0;
  for (let i = 0; i < arr.length; i++) {
    const rec = arr[i] || {};
    const sym = normalize_symbol(rec.sym || rec.label || rec.symbol);
    const z = Number.isFinite(rec.Z) ? (rec.Z | 0) : parse_int_safe(rec.Z, 0);
    if (sym === "H" || z === 1) n += 1;
  }
  return n;
}


function replacement_collect_coords_from_ocl(mol, atomLimit) {
  const coords = [];
  if (!mol || typeof mol.getAllAtoms !== "function") return coords;
  let count = parse_int_safe(mol.getAllAtoms(), 0);
  if (Number.isFinite(atomLimit) && atomLimit >= 0) count = Math.min(count, atomLimit | 0);
  for (let i = 0; i < count; i++) {
    let x = 0;
    let y = 0;
    let z = 0;
    try { if (typeof mol.getAtomX === "function") x = Number(mol.getAtomX(i)) || 0; } catch (_) {}
    try { if (typeof mol.getAtomY === "function") y = Number(mol.getAtomY(i)) || 0; } catch (_) {}
    try { if (typeof mol.getAtomZ === "function") z = Number(mol.getAtomZ(i)) || 0; } catch (_) {}
    coords.push({ x, y, z });
  }
  return coords;
}

function replacement_coords_have_3d_signal(coords, eps = 1e-4) {
  const arr = Array.isArray(coords) ? coords : [];
  if (!arr.length) return false;
  let minz = Number(arr[0] && arr[0].z) || 0;
  let maxz = minz;
  for (let i = 1; i < arr.length; i++) {
    const z = Number(arr[i] && arr[i].z) || 0;
    if (z < minz) minz = z;
    if (z > maxz) maxz = z;
  }
  return maxz - minz > eps;
}

function replacement_clone_ocl_molecule(mol) {
  if (!mol) return null;
  try {
    if (typeof mol.getCompactCopy === "function") return mol.getCompactCopy();
  } catch (_) {}
  try {
    if (typeof mol.toMolfile === "function") {
      const OCL = replacement_get_module();
      if (OCL && OCL.Molecule && typeof OCL.Molecule.fromMolfile === "function") {
        return OCL.Molecule.fromMolfile(String(mol.toMolfile() || ""));
      }
    }
  } catch (_) {}
  return null;
}

function replacement_collect_implicit_hydrogen_info(mol, atomLimit) {
  if (!mol || typeof mol.getImplicitHydrogens !== "function" || typeof mol.getAllAtoms !== "function") {
    return { supported: false, implicitHydrogenCount: 0, perAtom: [] };
  }
  let count = parse_int_safe(mol.getAllAtoms(), 0);
  if (Number.isFinite(atomLimit) && atomLimit >= 0) count = Math.min(count, atomLimit | 0);
  const perAtom = new Array(count);
  let implicitHydrogenCount = 0;
  for (let i = 0; i < count; i++) {
    try {
      const n = normalize_count_safe(mol.getImplicitHydrogens(i));
      if (n == null) return { supported: false, implicitHydrogenCount: 0, perAtom: [] };
      perAtom[i] = n;
      implicitHydrogenCount += n;
    } catch (_) {
      return { supported: false, implicitHydrogenCount: 0, perAtom: [] };
    }
  }
  return { supported: true, implicitHydrogenCount, perAtom };
}

function replacement_try_conformer_geometry(baseMol, atomLimit, warnings, diagnostics) {
  const OCL = replacement_get_module();
  if (!OCL || typeof OCL.ConformerGenerator !== "function") return null;
  const work = replacement_clone_ocl_molecule(baseMol);
  if (!work) {
    push_unique_warning(warnings, "replacement_conformer_clone_failed");
    return null;
  }

  try {
    const generator = new OCL.ConformerGenerator(1);
    const generated = generator.getOneConformerAsMolecule(work);
    const geomMol = generated || work;
    let geometryOrigin = "replacement_3d_conformer";

    if (OCL.ForceFieldMMFF94 && typeof OCL.ForceFieldMMFF94 === "function") {
      try {
        const tableName = OCL.ForceFieldMMFF94.MMFF94 || "MMFF94";
        const ff = new OCL.ForceFieldMMFF94(geomMol, tableName);
        const rc = typeof ff.minimise === "function" ? ff.minimise() : null;
        if (rc === 0) geometryOrigin = "replacement_3d_mmff94";
      } catch (_) {
        push_unique_warning(warnings, "replacement_mmff94_failed");
      }
    }

    const coords = replacement_collect_coords_from_ocl(geomMol, atomLimit);
    if (coords.length === atomLimit && !are_atoms_degenerate(coords) && replacement_coords_have_3d_signal(coords)) {
      if (diagnostics && typeof diagnostics === "object") diagnostics.strongGeometry = geometryOrigin;
      return { coords, geometryOrigin };
    }

    push_unique_warning(warnings, "replacement_conformer_failed");
    return null;
  } catch (_) {
    push_unique_warning(warnings, "replacement_conformer_failed");
    return null;
  }
}

function replacement_try_2d_geometry(baseMol, atomLimit, warnings) {
  const work = replacement_clone_ocl_molecule(baseMol);
  if (!work) return null;
  try {
    let coords = replacement_collect_coords_from_ocl(work, atomLimit);
    if (coords.length === atomLimit && !are_atoms_degenerate(coords)) {
      return { coords, geometryOrigin: replacement_coords_have_3d_signal(coords) ? "replacement_3d_existing" : "replacement_2d" };
    }
    if (typeof work.inventCoordinates === "function") work.inventCoordinates();
    coords = replacement_collect_coords_from_ocl(work, atomLimit);
    if (coords.length === atomLimit && !are_atoms_degenerate(coords)) {
      return { coords, geometryOrigin: replacement_coords_have_3d_signal(coords) ? "replacement_3d_existing" : "replacement_2d" };
    }
  } catch (_) {
    push_unique_warning(warnings, "replacement_inventCoordinates_failed");
  }
  return null;
}

function replacement_collect_atoms_from_ocl(mol) {
  const atoms = [];
  if (!mol || typeof mol.getAllAtoms !== "function") return atoms;
  const count = parse_int_safe(mol.getAllAtoms(), 0);
  for (let i = 0; i < count; i++) {
    let sym = "";
    let Z = 0;
    let x = 0;
    let y = 0;
    let z = 0;
    try { if (typeof mol.getAtomLabel === "function") sym = normalize_symbol(mol.getAtomLabel(i)); } catch (_) {}
    try { if (typeof mol.getAtomicNo === "function") Z = parse_int_safe(mol.getAtomicNo(i), 0); } catch (_) {}
    try { if (typeof mol.getAtomX === "function") x = Number(mol.getAtomX(i)) || 0; } catch (_) {}
    try { if (typeof mol.getAtomY === "function") y = Number(mol.getAtomY(i)) || 0; } catch (_) {}
    try { if (typeof mol.getAtomZ === "function") z = Number(mol.getAtomZ(i)) || 0; } catch (_) {}
    if (!(Z > 0) && sym) Z = local_symbol_to_Z(sym);
    if (!sym && Z > 0 && ELEMENT_SYMBOLS[Z]) sym = ELEMENT_SYMBOLS[Z];
    atoms.push({
      atom_idx: i,
      source_atom_id: i + 1,
      sym: sym || "",
      Z: Z > 0 ? Z : 0,
      x,
      y,
      z
    });
  }
  return atoms;
}

function replacement_collect_bonds_from_ocl(mol) {
  const bonds = [];
  if (!mol || typeof mol.getAllBonds !== "function") return bonds;
  const count = parse_int_safe(mol.getAllBonds(), 0);
  for (let i = 0; i < count; i++) {
    let a0 = 0;
    let a1 = 0;
    let order = 1;
    try { if (typeof mol.getBondAtom === "function") a0 = parse_int_safe(mol.getBondAtom(0, i), 0); } catch (_) {}
    try { if (typeof mol.getBondAtom === "function") a1 = parse_int_safe(mol.getBondAtom(1, i), 0); } catch (_) {}
    try { if (typeof mol.getBondOrder === "function") order = parse_bond_order_token(mol.getBondOrder(i)); } catch (_) {}
    bonds.push([a0, a1, order]);
  }
  return bonds;
}


function replacement_build_identity_from_atoms(atoms, meta) {
  const arr = Array.isArray(atoms) ? atoms : [];
  const zs = new Array(arr.length);
  const syms = new Array(arr.length);
  const uniq = new Set();
  for (let i = 0; i < arr.length; i++) {
    const atom = arr[i] || {};
    const sym = normalize_symbol(atom.sym);
    const z = Number.isFinite(atom.Z) ? (atom.Z | 0) : parse_int_safe(atom.Z, 0);
    zs[i] = z > 0 ? z : 0;
    syms[i] = sym || null;
    if (sym) uniq.add(sym);
  }

  const explicitHydrogenCount = replacement_count_h_atoms(arr);
  const implicitHydrogenSupported = !!(meta && meta.implicitHydrogenSupported);
  const implicitHydrogenCount = implicitHydrogenSupported && Number.isFinite(meta && meta.implicitHydrogenCount)
    ? Math.max(0, meta.implicitHydrogenCount | 0)
    : 0;

  let hydrogenMode = "unsupported";
  if (explicitHydrogenCount > 0 && !(implicitHydrogenSupported && implicitHydrogenCount > 0)) hydrogenMode = "explicit";
  else if (explicitHydrogenCount > 0 && implicitHydrogenSupported && implicitHydrogenCount > 0) hydrogenMode = "mixed";
  else if (implicitHydrogenSupported && implicitHydrogenCount > 0) hydrogenMode = "implicit_only";
  else if (implicitHydrogenSupported && explicitHydrogenCount === 0 && implicitHydrogenCount === 0) hydrogenMode = "none";

  const uniqueSymbols = Array.from(uniq).sort();
  const uiUniqueSymbols = hydrogenMode === "implicit_only" || hydrogenMode === "mixed"
    ? merge_unique_symbols(uniqueSymbols, "H")
    : uniqueSymbols.slice(0);

  return {
    Zs: zs,
    syms: syms,
    atomCount: arr.length,
    hydrogenCount: hydrogenMode === "implicit_only" || hydrogenMode === "mixed" ? (explicitHydrogenCount + implicitHydrogenCount) : explicitHydrogenCount,
    uniqueSymbols,
    uiUniqueSymbols,
    explicitHydrogenCount,
    implicitHydrogenCount,
    hydrogenMode,
    hasHydrogenElement: explicitHydrogenCount > 0 || (implicitHydrogenSupported && implicitHydrogenCount > 0),
    explicitHWorked: explicitHydrogenCount > 0 && !!(meta && meta.explicitHWorked) && !(implicitHydrogenSupported && implicitHydrogenCount > 0),
    source: (meta && meta.identitySource) || "replacement_backend",
    hydrogenInfoSource: implicitHydrogenSupported || explicitHydrogenCount > 0 ? ((meta && meta.hydrogenInfoSource) || "replacement_backend") : "unsupported",
    hydrogenInfoKind:
      hydrogenMode === "explicit"
        ? "explicit"
        : hydrogenMode === "mixed"
          ? "mixed"
          : hydrogenMode === "implicit_only"
            ? "implicit"
            : "none"
  };
}



function replacement_make_heavy_snapshot(atoms, bonds) {
  const atomArr = Array.isArray(atoms) ? atoms : [];
  const bondArr = Array.isArray(bonds) ? bonds : [];
  const heavyMap = new Map();
  const heavyAtoms = [];

  for (let i = 0; i < atomArr.length; i++) {
    const atom = atomArr[i] || {};
    const sym = normalize_symbol(atom.sym);
    const z = Number.isFinite(atom.Z) ? (atom.Z | 0) : parse_int_safe(atom.Z, 0);
    if (sym === "H" || z === 1) continue;
    const rec = { idx: i, sym: sym || (z > 0 && ELEMENT_SYMBOLS[z] ? ELEMENT_SYMBOLS[z] : "X"), z: z > 0 ? z : 0 };
    heavyMap.set(i, heavyAtoms.length);
    heavyAtoms.push(rec);
  }

  const heavySymbols = heavyAtoms.map((a) => a.sym).sort();
  const localNeighbors = new Array(heavyAtoms.length);
  for (let i = 0; i < localNeighbors.length; i++) localNeighbors[i] = [];
  const bondSignatures = [];

  for (let i = 0; i < bondArr.length; i++) {
    const bond = bondArr[i];
    if (!Array.isArray(bond) || bond.length < 2) continue;
    const a0 = parse_int_safe(bond[0], -1);
    const a1 = parse_int_safe(bond[1], -1);
    if (!heavyMap.has(a0) || !heavyMap.has(a1)) continue;
    const h0 = heavyAtoms[heavyMap.get(a0)];
    const h1 = heavyAtoms[heavyMap.get(a1)];
    const order = parse_bond_order_token(bond[2]);
    const pair = [h0.sym || "X", h1.sym || "X"].sort();
    bondSignatures.push(pair[0] + ":" + order + ":" + pair[1]);
    localNeighbors[heavyMap.get(a0)].push((h1.sym || "X") + ":" + order);
    localNeighbors[heavyMap.get(a1)].push((h0.sym || "X") + ":" + order);
  }

  const localSignatures = heavyAtoms.map((atom, idx) => {
    const neigh = Array.isArray(localNeighbors[idx]) ? localNeighbors[idx].slice(0).sort() : [];
    return (atom.sym || "X") + "|" + neigh.join(",");
  }).sort();

  return {
    heavyAtomCount: heavyAtoms.length,
    heavySymbols,
    heavyBondSignatures: bondSignatures.sort(),
    heavyLocalSignatures: localSignatures
  };
}

function replacement_same_sorted_arrays(a, b) {
  const aa = Array.isArray(a) ? a : [];
  const bb = Array.isArray(b) ? b : [];
  if (aa.length !== bb.length) return false;
  for (let i = 0; i < aa.length; i++) {
    if (String(aa[i]) !== String(bb[i])) return false;
  }
  return true;
}

function replacement_compare_heavy_snapshots(before, after) {
  const a = before && typeof before === "object" ? before : replacement_make_heavy_snapshot([], []);
  const b = after && typeof after === "object" ? after : replacement_make_heavy_snapshot([], []);
  if ((a.heavyAtomCount | 0) !== (b.heavyAtomCount | 0)) {
    return { ok: false, code: "replacement_heavy_atom_count_mismatch" };
  }
  if (!replacement_same_sorted_arrays(a.heavySymbols, b.heavySymbols)) {
    return { ok: false, code: "replacement_heavy_atom_symbol_mismatch" };
  }
  if (!replacement_same_sorted_arrays(a.heavyBondSignatures, b.heavyBondSignatures)) {
    return { ok: false, code: "replacement_heavy_bond_graph_mismatch" };
  }
  if (!replacement_same_sorted_arrays(a.heavyLocalSignatures, b.heavyLocalSignatures)) {
    return { ok: false, code: "replacement_heavy_local_signature_mismatch" };
  }
  return { ok: true, code: "ok" };
}

function replacement_geometry_score(geom) {
  const g = geom && typeof geom === "object" ? geom : null;
  const coords = g && Array.isArray(g.coords) ? g.coords : [];
  if (!coords.length) return 0;
  if (are_atoms_degenerate(coords)) return 1;
  if (replacement_coords_have_3d_signal(coords)) return 4;
  return 3;
}

function replacement_choose_geometry_candidate(current, candidate, warnings, diagnostics) {
  const cur = current && typeof current === "object" ? current : null;
  const next = candidate && typeof candidate === "object" ? candidate : null;
  if (!next || !Array.isArray(next.coords) || !next.coords.length) return cur;
  if (!cur || !Array.isArray(cur.coords) || !cur.coords.length) return next;
  const curScore = replacement_geometry_score(cur);
  const nextScore = replacement_geometry_score(next);
  if (nextScore > curScore) return next;
  if (nextScore === curScore) {
    const curCount = Array.isArray(cur.coords) ? cur.coords.length : 0;
    const nextCount = Array.isArray(next.coords) ? next.coords.length : 0;
    if (nextCount > curCount) return next;
  }
  if (nextScore < curScore) {
    push_unique_warning(warnings, "replacement_geometry_regressed_after_materialization");
    if (diagnostics && typeof diagnostics === "object") {
      diagnostics.geometryRegression = {
        kept: cur.geometryOrigin || "none",
        rejected: next.geometryOrigin || "none",
        keptScore: curScore,
        rejectedScore: nextScore
      };
    }
  }
  return cur;
}


function replacement_clear_non_fatal_materialization_warnings(warnings) {
  if (!Array.isArray(warnings) || !warnings.length) return;
  const transient = new Set([
    "replacement_explicit_h_materialization_unavailable",
    "replacement_explicit_h_materialization_incomplete",
    "replacement_explicit_h_bridge_failed",
    "replacement_explicit_h_roundtrip_failed",
    "replacement_heavy_atom_mismatch_after_bridge",
    "replacement_heavy_atom_mismatch_after_materialization",
    "replacement_heavy_atom_count_mismatch",
    "replacement_heavy_atom_symbol_mismatch",
    "replacement_heavy_bond_graph_mismatch",
    "replacement_heavy_local_signature_mismatch",
    "replacement_geometry_regressed_after_materialization"
  ]);
  for (let i = warnings.length - 1; i >= 0; i--) {
    if (transient.has(String(warnings[i] || "").trim())) warnings.splice(i, 1);
  }
}

function replacement_try_materialize_explicit_h_with_ocl(baseMol, baseSnapshot, warnings, diagnostics) {
  const work = replacement_clone_ocl_molecule(baseMol);
  if (!work) {
    push_unique_warning(warnings, "replacement_explicit_h_materialization_unavailable");
    return null;
  }
  if (typeof work.addImplicitHydrogens !== "function") {
    push_unique_warning(warnings, "replacement_explicit_h_materialization_unavailable");
    return null;
  }

  try {
    const preAtomCount = parse_int_safe(work.getAllAtoms && work.getAllAtoms(), 0);
    const preCoords = replacement_collect_coords_from_ocl(work, preAtomCount);
    if ((preCoords.length !== preAtomCount || are_atoms_degenerate(preCoords)) && typeof work.inventCoordinates === "function") {
      try { work.inventCoordinates(); } catch (_) {}
    }
    work.addImplicitHydrogens();
    const atoms = replacement_collect_atoms_from_ocl(work);
    const bonds = replacement_collect_bonds_from_ocl(work);
    const explicitHydrogenCount = replacement_count_h_atoms(atoms);
    const implicitHydrogenInfo = replacement_collect_implicit_hydrogen_info(work, atoms.length);
    const heavyCheck = replacement_compare_heavy_snapshots(baseSnapshot, replacement_make_heavy_snapshot(atoms, bonds));
    if (!heavyCheck.ok) {
      push_unique_warning(warnings, "replacement_heavy_atom_mismatch_after_materialization");
      push_unique_warning(warnings, heavyCheck.code);
      return null;
    }
    if (!(explicitHydrogenCount > 0)) {
      push_unique_warning(warnings, "replacement_explicit_h_materialization_unavailable");
      return null;
    }
    if (implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0) {
      push_unique_warning(warnings, "replacement_explicit_h_materialization_incomplete");
      return null;
    }

    let coords = replacement_collect_coords_from_ocl(work, atoms.length).map(clone_coord_record);
    const needRebuildCoords = coords.length !== atoms.length || !coords.length || are_atoms_degenerate(coords);
    if (needRebuildCoords && typeof work.inventCoordinates === "function") {
      try { work.inventCoordinates(); } catch (_) {}
      coords = replacement_collect_coords_from_ocl(work, atoms.length).map(clone_coord_record);
    }

    const geometryOrigin = coords.length === atoms.length && coords.length > 0 && !are_atoms_degenerate(coords)
      ? (replacement_coords_have_3d_signal(coords) ? "replacement_3d_existing" : "replacement_2d")
      : "replacement_none";

    dlog("replacement OCL materialized geometry", {
      atomCount: atoms.length,
      coordCount: coords.length,
      geometryOrigin,
      degenerate: are_atoms_degenerate(coords)
    });

    if (diagnostics && typeof diagnostics === "object") diagnostics.explicitMaterialization = "openchemlib";
    return {
      atoms,
      bonds,
      coords,
      geometryOrigin,
      explicitHydrogenCount,
      implicitHydrogenCount: 0,
      hydrogenAddedCount: Math.max(0, explicitHydrogenCount - replacement_count_h_atoms(replacement_collect_atoms_from_ocl(baseMol))),
      expandWorked: true,
      authoritativeTopologySource: "openchemlib",
      authoritativeHydrogenSource: "openchemlib",
      authoritativeGeometrySource: geometryOrigin === "replacement_none" ? "fallback" : "openchemlib",
      authoritativeMolKind: "openchemlib_materialized",
      identity: replacement_build_identity_from_atoms(atoms, {
        implicitHydrogenSupported: true,
        implicitHydrogenCount: 0,
        explicitHWorked: true,
        identitySource: "replacement_backend",
        hydrogenInfoSource: "openchemlib"
      })
    };
  } catch (_) {
    push_unique_warning(warnings, "replacement_explicit_h_materialization_unavailable");
    return null;
  }
}

function replacement_try_rdkit_mol_from_authoritative_representation(baseMol, diagnostics) {
  const representations = [];
  try {
    if (baseMol && typeof baseMol.toMolfileV3 === "function") {
      const text = String(baseMol.toMolfileV3() || "");
      if (text.trim()) representations.push({ kind: "molfile_v3", text });
    }
  } catch (_) {}
  try {
    if (baseMol && typeof baseMol.toMolfile === "function") {
      const text = String(baseMol.toMolfile() || "");
      if (text.trim()) representations.push({ kind: "molfile_v2000", text });
    }
  } catch (_) {}
  try {
    if (baseMol && typeof baseMol.toIsomericSmiles === "function") {
      const text = String(baseMol.toIsomericSmiles() || "");
      if (text.trim()) representations.push({ kind: "isomeric_smiles", text });
    }
  } catch (_) {}
  try {
    if (baseMol && typeof baseMol.toSmiles === "function") {
      const text = String(baseMol.toSmiles() || "");
      if (text.trim()) representations.push({ kind: "smiles", text });
    }
  } catch (_) {}

  for (let i = 0; i < representations.length; i++) {
    const rep = representations[i];
    const mol = raw_get_mol(rep.text);
    if (mol) {
      if (diagnostics && typeof diagnostics === "object") diagnostics.bridgeRepresentation = rep.kind;
      return { mol, kind: rep.kind };
    }
  }
  return null;
}

function replacement_try_materialize_explicit_h_with_rdkit_bridge(baseMol, baseSnapshot, warnings, diagnostics) {
  const bridgeInit = replacement_try_rdkit_mol_from_authoritative_representation(baseMol, diagnostics);
  if (!bridgeInit || !bridgeInit.mol) {
    push_unique_warning(warnings, "replacement_explicit_h_bridge_failed");
    return null;
  }

  let bridgeMol = bridgeInit.mol;
  let baseBridgeMol = bridgeMol;
  try {
    const hsRes = add_hs_safe(bridgeMol);
    if (!hsRes || !hsRes.apiAvailable) {
      push_unique_warning(warnings, "replacement_explicit_h_bridge_failed");
      return null;
    }
    if (hsRes.mol && hsRes.mol !== bridgeMol) bridgeMol = hsRes.mol;

    let atomsBonds = get_atoms_bonds_safe(bridgeMol);
    const atoms = Array.isArray(atomsBonds.atoms) ? atomsBonds.atoms.map(clone_atom_record) : [];
    const bonds = Array.isArray(atomsBonds.bonds) ? atomsBonds.bonds.map(clone_bond_record).filter(Boolean) : [];
    const explicitHydrogenCount = count_explicit_h_atoms_in_atoms(atoms);
    if (!(explicitHydrogenCount > 0)) {
      push_unique_warning(warnings, "replacement_explicit_h_bridge_failed");
      return null;
    }

    const heavyCheck = replacement_compare_heavy_snapshots(baseSnapshot, replacement_make_heavy_snapshot(atoms, bonds));
    if (!heavyCheck.ok) {
      push_unique_warning(warnings, "replacement_explicit_h_roundtrip_failed");
      push_unique_warning(warnings, "replacement_heavy_atom_mismatch_after_bridge");
      push_unique_warning(warnings, heavyCheck.code);
      return null;
    }

    let coords = get_molblock_coordinates_only(bridgeMol).map(clone_coord_record);
    let geometryOrigin = "replacement_none";
    const coordsUsable = coords.length === atoms.length && coords.length > 0 && !are_atoms_degenerate(coords);
    if (coordsUsable) {
      geometryOrigin = replacement_coords_have_3d_signal(coords) ? "replacement_3d_existing" : "replacement_2d";
    }

    if (!(coords.length === atoms.length && coords.length > 0 && replacement_coords_have_3d_signal(coords) && !are_atoms_degenerate(coords))) {
      const embedded = raw_embed(bridgeMol);
      if (embedded) raw_uff_optimize(bridgeMol);
      if (embedded) {
        const embedIdentity = replacement_build_identity_from_atoms(atoms, {
          implicitHydrogenSupported: true,
          implicitHydrogenCount: 0,
          explicitHWorked: true,
          identitySource: "replacement_backend",
          hydrogenInfoSource: "openchemlib_to_rdkit_bridge"
        });
        const snapshot = refresh_geometry_snapshot(bridgeMol, embedIdentity, bonds);
        const embedCoords = Array.isArray(snapshot.coords) ? snapshot.coords.map(clone_coord_record) : [];
        if (embedCoords.length === atoms.length && !are_atoms_degenerate(embedCoords) && replacement_coords_have_3d_signal(embedCoords)) {
          coords = embedCoords;
          geometryOrigin = "replacement_3d_mmff94";
        }
      }
    }

    if (!(coords.length === atoms.length && coords.length > 0) && (raw_compute2d(bridgeMol) || raw_generate_aligned_coords(bridgeMol))) {
      const identity2d = replacement_build_identity_from_atoms(atoms, {
        implicitHydrogenSupported: true,
        implicitHydrogenCount: 0,
        explicitHWorked: true,
        identitySource: "replacement_backend",
        hydrogenInfoSource: "openchemlib_to_rdkit_bridge"
      });
      const snapshot2d = refresh_geometry_snapshot(bridgeMol, identity2d, bonds);
      const coords2d = Array.isArray(snapshot2d.coords) ? snapshot2d.coords.map(clone_coord_record) : [];
      if (coords2d.length === atoms.length && coords2d.length > 0 && !are_atoms_degenerate(coords2d)) {
        coords = coords2d;
        geometryOrigin = replacement_coords_have_3d_signal(coords2d) ? "replacement_3d_existing" : "replacement_2d";
      }
    }

    if (diagnostics && typeof diagnostics === "object") diagnostics.explicitMaterialization = "openchemlib_to_rdkit_bridge";
    return {
      atoms,
      bonds,
      coords,
      geometryOrigin,
      explicitHydrogenCount,
      implicitHydrogenCount: 0,
      hydrogenAddedCount: Math.max(0, explicitHydrogenCount - replacement_count_h_atoms(replacement_collect_atoms_from_ocl(baseMol))),
      expandWorked: true,
      authoritativeTopologySource: "openchemlib",
      authoritativeHydrogenSource: "openchemlib_to_rdkit_bridge",
      authoritativeGeometrySource: geometryOrigin === "replacement_none" ? "fallback" : "rdkit",
      authoritativeMolKind: "rdkit_bridge_mol",
      bridgeRepresentationKind: bridgeInit.kind,
      identity: replacement_build_identity_from_atoms(atoms, {
        implicitHydrogenSupported: true,
        implicitHydrogenCount: 0,
        explicitHWorked: true,
        identitySource: "replacement_backend",
        hydrogenInfoSource: "openchemlib_to_rdkit_bridge"
      })
    };
  } finally {
    if (bridgeMol && bridgeMol !== baseBridgeMol) safe_delete(baseBridgeMol);
    safe_delete(bridgeMol);
  }
}


function build_rdkit_explicit_h_fallback_payload(baseMol) {
  if (!baseMol) return null;

  let workMol = baseMol;
  const before = get_atoms_bonds_safe(baseMol);
  const beforeAtoms = Array.isArray(before.atoms) ? before.atoms.map(clone_atom_record) : [];
  const beforeBonds = Array.isArray(before.bonds) ? before.bonds.map(clone_bond_record).filter(Boolean) : [];
  const baseSnapshot = replacement_make_heavy_snapshot(beforeAtoms, beforeBonds);

  const hsRes = add_hs_safe(workMol);
  if (!hsRes || !hsRes.apiAvailable || !hsRes.mol) return null;
  if (hsRes.mol && hsRes.mol !== workMol) workMol = hsRes.mol;

  let identity = replacement_build_identity_from_atoms(get_atoms_bonds_safe(workMol).atoms || [], {
    implicitHydrogenSupported: true,
    implicitHydrogenCount: 0,
    explicitHWorked: true,
    identitySource: "rdkit_backend",
    hydrogenInfoSource: "rdkit_addHs"
  });
  let snapshot = refresh_geometry_snapshot(workMol, identity, beforeBonds);
  let atoms = Array.isArray(snapshot.atoms) ? snapshot.atoms.map(clone_atom_record) : [];
  let bonds = Array.isArray(snapshot.bonds) ? snapshot.bonds.map(clone_bond_record).filter(Boolean) : [];
  let coords = Array.isArray(snapshot.coords) ? snapshot.coords.map(clone_coord_record) : [];

  const explicitHydrogenCount = count_explicit_h_atoms_in_atoms(atoms);
  if (!(explicitHydrogenCount > 0)) {
    if (workMol !== baseMol) safe_delete(workMol);
    return null;
  }

  const heavyCheck = replacement_compare_heavy_snapshots(baseSnapshot, replacement_make_heavy_snapshot(atoms, bonds));
  if (!heavyCheck.ok) {
    if (workMol !== baseMol) safe_delete(workMol);
    return null;
  }

  let coordsUsable = coords.length === atoms.length && coords.length > 0 && !are_atoms_degenerate(coords);
  let geometryOrigin = "rdkit";

  if (!coordsUsable) {
    const did2d = raw_compute2d(workMol) || raw_generate_aligned_coords(workMol);
    if (did2d) {
      identity = replacement_build_identity_from_atoms(atoms, {
        implicitHydrogenSupported: true,
        implicitHydrogenCount: 0,
        explicitHWorked: true,
        identitySource: "rdkit_backend",
        hydrogenInfoSource: "rdkit_addHs"
      });
      snapshot = refresh_geometry_snapshot(workMol, identity, bonds);
      atoms = Array.isArray(snapshot.atoms) ? snapshot.atoms.map(clone_atom_record) : atoms;
      bonds = Array.isArray(snapshot.bonds) ? snapshot.bonds.map(clone_bond_record).filter(Boolean) : bonds;
      coords = Array.isArray(snapshot.coords) ? snapshot.coords.map(clone_coord_record) : coords;
      coordsUsable = coords.length === atoms.length && coords.length > 0 && !are_atoms_degenerate(coords);
      if (coordsUsable) geometryOrigin = "2d_generated";
    }
  }

  const payload = {
    ok: atoms.length > 0,
    backendName: "rdkit",
    source: "rdkit_explicit_h_fallback",
    geometryOrigin,
    authoritativeTopologySource: "rdkit",
    authoritativeHydrogenSource: "rdkit_addHs",
    authoritativeGeometrySource: coordsUsable ? "rdkit" : "fallback",
    authoritativeMolKind: "rdkit_explicit_h_mol",
    explicitHydrogenCount,
    implicitHydrogenCount: 0,
    baseAtomCount: beforeAtoms.length,
    baseBondCount: beforeBonds.length,
    hydrogenAddedCount: Math.max(0, explicitHydrogenCount - count_explicit_h_atoms_in_atoms(beforeAtoms)),
    expandWorked: !!hsRes.changed || explicitHydrogenCount > 0,
    usedWrapperMol: false,
    warnings: [],
    diagnostics: {
      explicitMaterialization: "rdkit_direct_addHs"
    },
    identity: replacement_build_identity_from_atoms(atoms, {
      implicitHydrogenSupported: true,
      implicitHydrogenCount: 0,
      explicitHWorked: true,
      identitySource: "rdkit_backend",
      hydrogenInfoSource: "rdkit_addHs"
    }),
    atoms,
    bonds,
    coords
  };

  return { mol: workMol, payload };
}


function replacement_extract_smiles_result(smiles) {
  const OCL = replacement_get_module();
  if (!OCL || !OCL.Molecule || typeof OCL.Molecule.fromSmiles !== "function") return null;

  let mol = null;
  const warnings = [];
  try {
    mol = OCL.Molecule.fromSmiles(String(smiles == null ? "" : smiles));
  } catch (e) {
    dlog("replacement backend parse failed", e);
    return null;
  }
  if (!mol) return null;

  let baseAtomCount = 0;
  let baseBondCount = 0;
  try { if (typeof mol.getAllAtoms === "function") baseAtomCount = parse_int_safe(mol.getAllAtoms(), 0); } catch (_) {}
  try { if (typeof mol.getAllBonds === "function") baseBondCount = parse_int_safe(mol.getAllBonds(), 0); } catch (_) {}

  let atoms = replacement_collect_atoms_from_ocl(mol);
  let bonds = replacement_collect_bonds_from_ocl(mol);
  let explicitHydrogenCount = replacement_count_h_atoms(atoms);
  let implicitHydrogenInfo = replacement_collect_implicit_hydrogen_info(mol, baseAtomCount);

  const diagnostics = {
    backendName: REPLACEMENT_BACKEND_NAME,
    replacementSource: replacement_get_status().source || REPLACEMENT_BACKEND_NAME,
    baseAtomCount,
    baseBondCount,
    explicitHydrogenCount,
    implicitHydrogenSupported: !!implicitHydrogenInfo.supported,
    implicitHydrogenCount: implicitHydrogenInfo.supported ? (implicitHydrogenInfo.implicitHydrogenCount | 0) : 0,
    strongGeometry: "none",
    authoritativeTopologySource: "openchemlib",
    authoritativeHydrogenSource: explicitHydrogenCount > 0 && !(implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0)
      ? "openchemlib"
      : implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0
        ? "implicit_only"
        : "fallback",
    authoritativeGeometrySource: "fallback",
    authoritativeMolKind: "openchemlib_molecule"
  };

  let geometry = null;
  const directCoords = replacement_collect_coords_from_ocl(mol, baseAtomCount);
  if (directCoords.length === baseAtomCount && !are_atoms_degenerate(directCoords)) {
    geometry = {
      coords: directCoords.map(clone_coord_record),
      geometryOrigin: replacement_coords_have_3d_signal(directCoords) ? "replacement_3d_existing" : "replacement_2d"
    };
  }

  if (!geometry || geometry.geometryOrigin === "replacement_2d") {
    const strongGeometry = replacement_try_conformer_geometry(mol, baseAtomCount, warnings, diagnostics);
    if (strongGeometry && strongGeometry.coords && strongGeometry.coords.length === baseAtomCount) {
      geometry = replacement_choose_geometry_candidate(geometry, strongGeometry, warnings, diagnostics);
    }
  }

  if (!geometry) {
    geometry = replacement_try_2d_geometry(mol, baseAtomCount, warnings);
  }

  const baseHeavySnapshot = replacement_make_heavy_snapshot(atoms, bonds);
  let chosenMaterialization = null;
  const alreadyExplicitOnly = explicitHydrogenCount > 0 && !(implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0);
  const implicitHydrogenKnown = !!(implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0);
  const knownZeroHydrogen = !!(implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount === 0 && explicitHydrogenCount === 0);
  const shouldProbeExplicitHydrogens = !alreadyExplicitOnly && !knownZeroHydrogen && (implicitHydrogenKnown || explicitHydrogenCount === 0);
  if (shouldProbeExplicitHydrogens) {
    diagnostics.explicitHydrogenProbeReason = implicitHydrogenKnown ? "implicit_h_reported" : "identity_not_explicit";
    const oclMaterialized = replacement_try_materialize_explicit_h_with_ocl(mol, baseHeavySnapshot, warnings, diagnostics);
    if (oclMaterialized) chosenMaterialization = oclMaterialized;

    const bridgeMaterialized = replacement_try_materialize_explicit_h_with_rdkit_bridge(mol, baseHeavySnapshot, warnings, diagnostics);
    if (bridgeMaterialized) {
      if (!chosenMaterialization) chosenMaterialization = bridgeMaterialized;
      else {
        const chosenGeom = { coords: chosenMaterialization.coords, geometryOrigin: chosenMaterialization.geometryOrigin };
        const betterGeom = replacement_choose_geometry_candidate(chosenGeom, { coords: bridgeMaterialized.coords, geometryOrigin: bridgeMaterialized.geometryOrigin }, warnings, diagnostics);
        if (betterGeom && betterGeom.geometryOrigin === bridgeMaterialized.geometryOrigin) {
          chosenMaterialization = bridgeMaterialized;
        }
      }
    }

    if (!chosenMaterialization && implicitHydrogenKnown) {
      push_unique_warning(warnings, "replacement_explicit_h_materialization_unavailable");
    }
  }

  let authoritativeTopologySource = diagnostics.authoritativeTopologySource || "openchemlib";
  let authoritativeHydrogenSource = diagnostics.authoritativeHydrogenSource || "fallback";
  let authoritativeGeometrySource = diagnostics.authoritativeGeometrySource || "fallback";
  let authoritativeMolKind = diagnostics.authoritativeMolKind || "openchemlib_molecule";
  let hydrogenAddedCount = 0;
  let expandWorked = false;

  if (chosenMaterialization) {
    if (warnings.length) diagnostics.materializationAttemptWarnings = warnings.slice(0);
    replacement_clear_non_fatal_materialization_warnings(warnings);
    atoms = Array.isArray(chosenMaterialization.atoms) ? chosenMaterialization.atoms.map(clone_atom_record) : atoms;
    bonds = Array.isArray(chosenMaterialization.bonds) ? chosenMaterialization.bonds.map(clone_bond_record).filter(Boolean) : bonds;
    explicitHydrogenCount = Number.isFinite(chosenMaterialization.explicitHydrogenCount) ? (chosenMaterialization.explicitHydrogenCount | 0) : explicitHydrogenCount;
    hydrogenAddedCount = Number.isFinite(chosenMaterialization.hydrogenAddedCount) ? (chosenMaterialization.hydrogenAddedCount | 0) : Math.max(0, explicitHydrogenCount);
    implicitHydrogenInfo = { supported: true, implicitHydrogenCount: 0, perAtom: [] };
    expandWorked = !!chosenMaterialization.expandWorked;
    authoritativeTopologySource = chosenMaterialization.authoritativeTopologySource || authoritativeTopologySource;
    authoritativeHydrogenSource = chosenMaterialization.authoritativeHydrogenSource || authoritativeHydrogenSource;
    authoritativeMolKind = chosenMaterialization.authoritativeMolKind || authoritativeMolKind;
    if (Array.isArray(chosenMaterialization.coords) && chosenMaterialization.coords.length) {
      geometry = replacement_choose_geometry_candidate(geometry, {
        coords: chosenMaterialization.coords.map(clone_coord_record),
        geometryOrigin: chosenMaterialization.geometryOrigin || "replacement_none"
      }, warnings, diagnostics);
      authoritativeGeometrySource = chosenMaterialization.authoritativeGeometrySource || authoritativeGeometrySource;
    }
  }

  const coords = geometry && Array.isArray(geometry.coords) ? geometry.coords.map(clone_coord_record) : [];
  const geometryOrigin = geometry ? (geometry.geometryOrigin || "replacement_none") : "replacement_none";

  if (geometryOrigin !== "replacement_none") {
    authoritativeGeometrySource = replacement_coords_have_3d_signal(coords)
      ? (geometryOrigin === "rdkit" ? "rdkit" : authoritativeGeometrySource === "fallback" ? "openchemlib" : authoritativeGeometrySource)
      : authoritativeGeometrySource === "fallback" ? "openchemlib" : authoritativeGeometrySource;
  }

  const identity = chosenMaterialization && chosenMaterialization.identity
    ? clone_identity_info(chosenMaterialization.identity)
    : replacement_build_identity_from_atoms(atoms, {
        implicitHydrogenSupported: !!implicitHydrogenInfo.supported,
        implicitHydrogenCount: implicitHydrogenInfo.supported ? (implicitHydrogenInfo.implicitHydrogenCount | 0) : 0,
        explicitHWorked: explicitHydrogenCount > 0 && !(implicitHydrogenInfo.supported && implicitHydrogenInfo.implicitHydrogenCount > 0),
        identitySource: "replacement_backend",
        hydrogenInfoSource: authoritativeHydrogenSource
      });

  diagnostics.explicitHydrogenCount = explicitHydrogenCount;
  diagnostics.implicitHydrogenCount = implicitHydrogenInfo.supported ? (implicitHydrogenInfo.implicitHydrogenCount | 0) : 0;
  diagnostics.authoritativeTopologySource = authoritativeTopologySource;
  diagnostics.authoritativeHydrogenSource = authoritativeHydrogenSource;
  diagnostics.authoritativeGeometrySource = authoritativeGeometrySource;
  diagnostics.authoritativeMolKind = authoritativeMolKind;
  diagnostics.hydrogenAddedCount = hydrogenAddedCount;

  return {
    ok: atoms.length > 0,
    backendName: REPLACEMENT_BACKEND_NAME,
    source: "replacement_backend",
    geometryOrigin,
    authoritativeTopologySource,
    authoritativeHydrogenSource,
    authoritativeGeometrySource,
    authoritativeMolKind,
    explicitHydrogenCount,
    implicitHydrogenCount: implicitHydrogenInfo.supported ? (implicitHydrogenInfo.implicitHydrogenCount | 0) : 0,
    baseAtomCount,
    baseBondCount,
    hydrogenAddedCount,
    expandWorked,
    warnings,
    diagnostics,
    identity,
    atoms,
    bonds,
    coords
  };
}


function req() {
  if (!_RDKit) throw new Error("RDKitModule ще не ініціалізований");
  return _RDKit;
}

function dev_enabled() {
  try {
    return !!(typeof window !== "undefined" && window && window.TEM_DEV === true);
  } catch (_) {
    return false;
  }
}

function dlog(...args) {
  if (dev_enabled()) console.log("[RDKit adapter]", ...args);
}

function normalize_symbol(sym) {
  let s = String(sym == null ? "" : sym).trim();
  if (!s) return "";
  if (s.length === 1) return s.toUpperCase();
  return s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
}

function local_symbol_to_Z(sym) {
  const norm = normalize_symbol(sym);
  return norm && SYMBOL_TO_Z[norm] ? SYMBOL_TO_Z[norm] : 0;
}

function safe_delete(obj) {
  try {
    if (obj && typeof obj.delete === "function") obj.delete();
  } catch (_) {}
}

function raw_method(obj, names) {
  if (!obj) return null;
  for (let i = 0; i < names.length; i++) {
    const fn = obj[names[i]];
    if (typeof fn === "function") return fn.bind(obj);
  }
  return null;
}

function parse_int_safe(v, fallback = 0) {
  const n = parseInt(v, 10);
  return Number.isFinite(n) ? n : fallback;
}

function parse_float_safe(v, fallback = NaN) {
  const n = parseFloat(v);
  return Number.isFinite(n) ? n : fallback;
}

function parse_v2000_count_segment(line, start, end) {
  const seg = String(line || "").slice(start, end);
  const n = parseInt(seg, 10);
  if (Number.isFinite(n)) return n;
  return null;
}

function parse_v2000_counts(lines) {
  const counts = String((lines && lines[3]) || "");
  let atomCount = parse_v2000_count_segment(counts, 0, 3);
  let bondCount = parse_v2000_count_segment(counts, 3, 6);

  if (!Number.isFinite(atomCount)) {
    const m = counts.match(/^\s*(\d+)\s+(\d+)/);
    atomCount = m ? parseInt(m[1], 10) : 0;
    bondCount = m ? parseInt(m[2], 10) : 0;
  }

  atomCount = Number.isFinite(atomCount) ? atomCount : 0;
  bondCount = Number.isFinite(bondCount) ? bondCount : 0;
  return { atomCount, bondCount };
}

function v2000_atom_symbol_from_line(line) {
  const ln = String(line || "");
  const fixed = normalize_symbol(ln.slice(31, 34));
  if (fixed) return fixed;
  const parts = ln.trim().split(/\s+/);
  return normalize_symbol(parts[3] || "");
}

function v2000_atom_coords_from_line(line) {
  const ln = String(line || "");
  let x = parse_float_safe(ln.slice(0, 10));
  let y = parse_float_safe(ln.slice(10, 20));
  let z = parse_float_safe(ln.slice(20, 30));

  if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
    const parts = ln.trim().split(/\s+/);
    x = parse_float_safe(parts[0]);
    y = parse_float_safe(parts[1]);
    z = parse_float_safe(parts[2]);
  }

  return {
    x: Number.isFinite(x) ? x : 0,
    y: Number.isFinite(y) ? y : 0,
    z: Number.isFinite(z) ? z : 0
  };
}

function parse_bond_order_token(token) {
  const raw = String(token == null ? "" : token).trim();
  if (!raw) return 1;
  const up = raw.toUpperCase();
  if (up === "AROMATIC" || up === "AROM") return 4;
  const n = parseFloat(raw);
  if (!Number.isFinite(n)) return 1;
  return n > 0 ? n : 1;
}

function parse_v2000_bond_line(line) {
  const ln = String(line || "");
  let i = parse_int_safe(ln.slice(0, 3), 0) - 1;
  let j = parse_int_safe(ln.slice(3, 6), 0) - 1;
  let order = parse_bond_order_token(ln.slice(6, 9));

  if (i < 0 || j < 0) {
    const parts = ln.trim().split(/\s+/);
    i = parse_int_safe(parts[0], 0) - 1;
    j = parse_int_safe(parts[1], 0) - 1;
    order = parse_bond_order_token(parts[2]);
  }

  if (i < 0 || j < 0) return null;
  return [i, j, order];
}

function collect_v3000_block(lines, beginToken, endToken) {
  const out = [];
  let inBlock = false;
  for (let i = 0; i < lines.length; i++) {
    const raw = String(lines[i] || "");
    const ln = raw.trim();
    if (!ln) continue;
    if (ln.startsWith(beginToken)) {
      inBlock = true;
      continue;
    }
    if (ln.startsWith(endToken)) break;
    if (!inBlock || !ln.startsWith("M  V30 ")) continue;
    out.push(ln.slice(7).trim());
  }
  return out;
}

function parse_v3000_atom_line(line, atomIdx) {
  const parts = String(line || "").trim().split(/\s+/);
  if (parts.length < 5) return null;
  const sym = normalize_symbol(parts[1]);
  const x = parse_float_safe(parts[2], 0);
  const y = parse_float_safe(parts[3], 0);
  const z = parse_float_safe(parts[4], 0);
  return {
    atom_idx: atomIdx,
    source_atom_id: parse_int_safe(parts[0], atomIdx + 1),
    sym,
    Z: 0,
    x,
    y,
    z
  };
}

function parse_v3000_bond_line(line, atomIdToIdx) {
  const parts = String(line || "").trim().split(/\s+/);
  if (parts.length < 4) return null;
  const order = parse_bond_order_token(parts[1]);
  const beginId = parse_int_safe(parts[2], 0);
  const endId = parse_int_safe(parts[3], 0);
  if (!beginId || !endId) return null;
  const i = atomIdToIdx.has(beginId) ? atomIdToIdx.get(beginId) : beginId - 1;
  const j = atomIdToIdx.has(endId) ? atomIdToIdx.get(endId) : endId - 1;
  if (i < 0 || j < 0) return null;
  return [i, j, order];
}

function summarize_parsed_molblock(parsed) {
  const seen = new Set();
  let hydrogenCount = 0;
  const atoms = Array.isArray(parsed && parsed.atoms) ? parsed.atoms : [];
  for (let i = 0; i < atoms.length; i++) {
    const sym = normalize_symbol(atoms[i] && atoms[i].sym);
    if (!sym) continue;
    seen.add(sym);
    if (sym === "H") hydrogenCount += 1;
  }
  return {
    atomCount: atoms.length,
    hydrogenCount,
    uniqueSymbols: Array.from(seen).sort(),
    format: parsed && parsed.format ? parsed.format : "none"
  };
}

function parse_molblock_text(molblock) {
  const text = String(molblock == null ? "" : molblock);
  if (!text.trim()) {
    return {
      atoms: [],
      bonds: [],
      format: "none"
    };
  }

  const lines = text.replace(/\r/g, "").split("\n");
  const hasV3000 = lines.some((ln, idx) => idx < 24 && String(ln || "").includes("V3000"));

  if (hasV3000) {
    const atomLines = collect_v3000_block(lines, "M  V30 BEGIN ATOM", "M  V30 END ATOM");
    const atoms = [];
    const atomIdToIdx = new Map();
    for (let i = 0; i < atomLines.length; i++) {
      const atom = parse_v3000_atom_line(atomLines[i], atoms.length);
      if (!atom) continue;
      atoms.push(atom);
      atomIdToIdx.set(atom.source_atom_id, atom.atom_idx);
    }

    const bondLines = collect_v3000_block(lines, "M  V30 BEGIN BOND", "M  V30 END BOND");
    const bonds = [];
    for (let i = 0; i < bondLines.length; i++) {
      const bond = parse_v3000_bond_line(bondLines[i], atomIdToIdx);
      if (bond) bonds.push(bond);
    }

    return {
      atoms,
      bonds,
      format: "V3000"
    };
  }

  const counts = parse_v2000_counts(lines);
  const atoms = [];
  for (let i = 0; i < counts.atomCount; i++) {
    const line = lines[4 + i] || "";
    const coords = v2000_atom_coords_from_line(line);
    atoms.push({
      atom_idx: i,
      source_atom_id: i + 1,
      sym: v2000_atom_symbol_from_line(line),
      Z: 0,
      x: coords.x,
      y: coords.y,
      z: coords.z
    });
  }

  const bonds = [];
  const bondStart = 4 + counts.atomCount;
  for (let k = 0; k < counts.bondCount; k++) {
    const bond = parse_v2000_bond_line(lines[bondStart + k] || "");
    if (bond) bonds.push(bond);
  }

  return {
    atoms,
    bonds,
    format: "V2000"
  };
}

function raw_ptable_getAtomicNumber(sym) {
  try {
    const mod = req();
    const norm = normalize_symbol(sym);
    if (!norm) return 0;
    if (typeof mod.getAtomicNumber === "function") {
      const z = mod.getAtomicNumber(norm);
      return Number.isFinite(z) ? (z | 0) : 0;
    }
  } catch (_) {}
  return 0;
}

function symbol_to_Z_safe(sym) {
  const norm = normalize_symbol(sym);
  if (!norm) return 0;
  const rdkitZ = raw_ptable_getAtomicNumber(norm);
  if (rdkitZ > 0) return rdkitZ;
  return local_symbol_to_Z(norm);
}

function raw_get_mol(smiles) {
  try {
    const mod = req();
    if (typeof mod.get_mol === "function") return mod.get_mol(smiles);
    if (mod.Mol) return new mod.Mol(smiles);
  } catch (_) {}
  return null;
}

function raw_get_molblock(mol) {
  try {
    const fn = raw_method(mol, [
      "get_molblock",
      "MolToMolBlock",
      "getMolblock",
      "getMolBlock"
    ]);
    return fn ? String(fn() || "") : "";
  } catch (_) {
    return "";
  }
}

function raw_get_json(mol) {
  try {
    const fn = raw_method(mol, [
      "get_json",
      "getJson",
      "MolToJSON",
      "molToJSON"
    ]);
    return fn ? String(fn() || "") : "";
  } catch (_) {
    return "";
  }
}

function raw_get_num_atoms(mol) {
  try {
    const fn = raw_method(mol, ["__tem_native_getNumAtoms", "getNumAtoms"]);
    if (!fn) return 0;
    const n = fn();
    return Number.isFinite(n) ? (n | 0) : 0;
  } catch (_) {
    return 0;
  }
}

function raw_get_num_bonds(mol) {
  try {
    const fn = raw_method(mol, ["__tem_native_getNumBonds", "getNumBonds"]);
    if (!fn) return 0;
    const n = fn();
    return Number.isFinite(n) ? (n | 0) : 0;
  } catch (_) {
    return 0;
  }
}

function raw_get_atom_with_idx(mol, idx) {
  try {
    const fn = raw_method(mol, ["get_atom_with_idx", "getAtomWithIdx"]);
    return fn ? fn(idx) : null;
  } catch (_) {
    return null;
  }
}

function raw_atom_get_atomic_number(atom) {
  try {
    const fn = raw_method(atom, ["get_atomic_num", "getAtomicNum"]);
    if (!fn) return 0;
    const z = fn();
    return Number.isFinite(z) ? (z | 0) : 0;
  } catch (_) {
    return 0;
  }
}

function raw_atom_get_symbol(atom) {
  try {
    const fn = raw_method(atom, ["get_symbol", "getSymbol"]);
    return fn ? normalize_symbol(fn()) : "";
  } catch (_) {
    return "";
  }
}

function raw_addHs(mol) {
  if (!mol) return null;
  try {
    const fn = raw_method(mol, [
      "__tem_native_addHs",
      "__tem_native_add_hs",
      "__tem_native_AddHs",
      "addHs",
      "AddHs"
    ]);
    return fn ? fn() : null;
  } catch (_) {
    return null;
  }
}

function raw_embed(mol) {
  if (!mol) return false;
  try {
    const fn = raw_method(mol, [
      "__tem_native_embedMolecule",
      "__tem_native_EmbedMolecule",
      "embedMolecule",
      "EmbedMolecule"
    ]);
    if (!fn) return false;
    const out = fn();
    if (typeof out === "boolean") return out;
    if (Number.isFinite(out)) return out === 0 || out === 1 ? !!out : true;
    return !!out;
  } catch (_) {
    return false;
  }
}

function raw_uff_optimize(mol) {
  if (!mol) return false;
  try {
    const fn = raw_method(mol, [
      "__tem_native_UFFOptimizeMolecule",
      "__tem_native_uffOptimizeMolecule",
      "UFFOptimizeMolecule",
      "uffOptimizeMolecule"
    ]);
    if (!fn) return false;
    fn();
    return true;
  } catch (_) {
    return false;
  }
}

function raw_compute2d(mol) {
  if (!mol) return false;
  try {
    const fn = raw_method(mol, ["compute2DCoords", "compute_2d_coords"]);
    if (!fn) return false;
    fn();
    return true;
  } catch (_) {
    return false;
  }
}

function raw_generate_aligned_coords(mol) {
  if (!mol) return false;
  try {
    const fn = raw_method(mol, [
      "__tem_native_generate_alignedCoords",
      "__tem_native_generateAlignedCoords",
      "generate_aligned_coords",
      "generate_alignedCoords",
      "generateAlignedCoords"
    ]);
    if (!fn) return false;
    fn();
    return true;
  } catch (_) {
    return false;
  }
}

function get_molblock_text(mol) {
  return raw_get_molblock(mol);
}

function get_json_text(mol) {
  return raw_get_json(mol);
}


function get_num_atoms_safe(mol) {
  const replacement = replacement_get_cache(mol);
  if (replacement && Array.isArray(replacement.atoms)) return replacement.atoms.length;
  const rawCount = raw_get_num_atoms(mol);
  if (rawCount > 0) return rawCount;
  return summarize_molblock(mol).atomCount;
}

function get_num_bonds_safe(mol) {
  const replacement = replacement_get_cache(mol);
  if (replacement && Array.isArray(replacement.bonds)) return replacement.bonds.length;
  const rawCount = raw_get_num_bonds(mol);
  if (rawCount > 0) return rawCount;
  const parsed = parse_molblock_text(get_molblock_text(mol));
  return Array.isArray(parsed.bonds) ? parsed.bonds.length : 0;
}

function summarize_molblock(mol) {
  const parsed = parse_molblock_text(get_molblock_text(mol));
  const summary = summarize_parsed_molblock(parsed);
  dlog("molblock summary", summary);
  return summary;
}

function decorate_atoms_with_Z(atoms) {
  for (let i = 0; i < atoms.length; i++) {
    const atom = atoms[i];
    atom.sym = normalize_symbol(atom.sym);
    atom.Z = symbol_to_Z_safe(atom.sym);
    atom.atom_idx = Number.isFinite(atom.atom_idx) ? atom.atom_idx : i;
  }
  return atoms;
}

function get_molblock_coordinates_only(mol) {
  const parsed = parse_molblock_text(get_molblock_text(mol));
  const parsedCoords = [];
  const parsedAtoms = Array.isArray(parsed.atoms) ? parsed.atoms : [];
  for (let i = 0; i < parsedAtoms.length; i++) {
    const a = parsedAtoms[i];
    if (!Number.isFinite(a.x) || !Number.isFinite(a.y) || !Number.isFinite(a.z)) continue;
    parsedCoords.push({ x: a.x, y: a.y, z: a.z });
  }
  return parsedCoords;
}

function get_coordinates_safe(mol) {
  const parsedCoords = get_molblock_coordinates_only(mol);

  const replacement = replacement_get_cache(mol);
  if (replacement) {
    if (parsedCoords.length > 0 && Array.isArray(replacement.atoms) && parsedCoords.length === replacement.atoms.length) {
      return parsedCoords;
    }
    if (Array.isArray(replacement.coords) && replacement.coords.length > 0) {
      return replacement.coords.map(clone_coord_record);
    }
  }

  return parsedCoords;
}

function get_atoms_bonds_safe(mol) {
  const replacement = replacement_get_cache(mol);
  if (replacement) {
    const atoms = Array.isArray(replacement.atoms) ? replacement.atoms.map(clone_atom_record) : [];
    const bonds = Array.isArray(replacement.bonds) ? replacement.bonds.map(clone_bond_record).filter(Boolean) : [];
    const coords = get_coordinates_safe(mol);
    if (coords.length === atoms.length) {
      for (let i = 0; i < atoms.length; i++) {
        atoms[i].x = coords[i].x;
        atoms[i].y = coords[i].y;
        atoms[i].z = coords[i].z;
      }
    }

    const uniqueSymbols = Array.from(new Set(atoms.map((a) => a.sym).filter(Boolean))).sort();
    dlog("atoms+bonds", {
      atomCount: atoms.length,
      bondCount: bonds.length,
      uniqueSymbols,
      source: replacement.source || "replacement_backend"
    });

    return { atoms, bonds };
  }

  const parsed = parse_molblock_text(get_molblock_text(mol));
  const atoms = decorate_atoms_with_Z(parsed.atoms || []);
  const bonds = Array.isArray(parsed.bonds) ? parsed.bonds : [];

  const uniqueSymbols = Array.from(new Set(atoms.map((a) => a.sym).filter(Boolean))).sort();
  dlog("atoms+bonds", {
    atomCount: atoms.length,
    bondCount: bonds.length,
    uniqueSymbols
  });

  return { atoms, bonds };
}

function add_hs_safe(mol) {
  const replacement = replacement_get_cache(mol);
  if (replacement) {
    const atomCountAfter = Array.isArray(replacement.atoms) ? replacement.atoms.length : 0;
    const explicitHydrogenCount = replacement.explicitHydrogenCount | 0;
    const atomCountBefore = Number.isFinite(replacement.baseAtomCount) ? (replacement.baseAtomCount | 0) : atomCountAfter;
    const result = {
      mol,
      changed: false,
      apiAvailable: false,
      returnedNewMol: false,
      atomCountBefore,
      atomCountAfter,
      molblockHydrogenCountAfter: explicitHydrogenCount,
      source: replacement.backendName || REPLACEMENT_BACKEND_NAME
    };
    dlog("AddHs result", result);
    return result;
  }

  const atomCountBefore = get_num_atoms_safe(mol);
  const apiFn = raw_method(mol, [
    "__tem_native_addHs",
    "__tem_native_add_hs",
    "__tem_native_AddHs",
    "addHs",
    "AddHs"
  ]);
  const apiAvailable = typeof apiFn === "function";

  let returned = null;
  if (apiAvailable) {
    try {
      returned = apiFn();
    } catch (_) {
      returned = null;
    }
  }

  let outMol = mol;
  let returnedNewMol = false;
  if (returned && returned !== mol) {
    outMol = returned;
    returnedNewMol = true;
  } else if (returned) {
    outMol = returned;
  }

  const atomCountAfter = get_num_atoms_safe(outMol);
  const mbSummaryAfter = summarize_molblock(outMol);
  const changed = atomCountAfter !== atomCountBefore || mbSummaryAfter.hydrogenCount > 0;

  const result = {
    mol: outMol,
    changed,
    apiAvailable,
    returnedNewMol,
    atomCountBefore,
    atomCountAfter,
    molblockHydrogenCountAfter: mbSummaryAfter.hydrogenCount
  };

  dlog("AddHs result", result);
  return result;
}

function safe_get_atom_atomicnum(mol, idx) {
  const atom = raw_get_atom_with_idx(mol, idx);
  if (atom) {
    const z = raw_atom_get_atomic_number(atom);
    if (z > 0) {
      safe_delete(atom);
      return z;
    }
    const sym = raw_atom_get_symbol(atom);
    safe_delete(atom);
    const safeZ = symbol_to_Z_safe(sym);
    if (safeZ > 0) return safeZ;
  }

  const extracted = get_atoms_bonds_safe(mol);
  const rec = extracted.atoms[idx];
  return rec && Number.isFinite(rec.Z) ? (rec.Z | 0) : 0;
}

function safe_get_bond_record(mol, idx) {
  const extracted = get_atoms_bonds_safe(mol);
  return extracted.bonds[idx] || null;
}


function normalize_count_safe(v) {
  const n = Number(v);
  if (!Number.isFinite(n)) return null;
  if (n < 0) return 0;
  return Math.max(0, Math.round(n));
}

function raw_atom_get_hcount(atom, names) {
  try {
    const fn = raw_method(atom, names);
    if (!fn) return null;
    return normalize_count_safe(fn());
  } catch (_) {
    return null;
  }
}

function collect_atom_api_hydrogen_info(mol, atomCount) {
  let explicitSum = 0;
  let implicitSum = 0;
  let totalSum = 0;
  let foundExplicit = false;
  let foundImplicit = false;
  let foundTotal = false;

  for (let i = 0; i < atomCount; i++) {
    const atom = raw_get_atom_with_idx(mol, i);
    if (!atom) continue;

    const exp = raw_atom_get_hcount(atom, ["getNumExplicitHs", "get_num_explicit_hs"]);
    const imp = raw_atom_get_hcount(atom, ["getNumImplicitHs", "get_num_implicit_hs"]);
    const total = raw_atom_get_hcount(atom, ["getTotalNumHs", "get_total_num_hs"]);

    if (exp != null) { explicitSum += exp; foundExplicit = true; }
    if (imp != null) { implicitSum += imp; foundImplicit = true; }
    if (total != null) { totalSum += total; foundTotal = true; }

    safe_delete(atom);
  }

  if (!foundExplicit && !foundImplicit && !foundTotal) {
    return { supported: false, explicitHydrogenCount: 0, implicitHydrogenCount: 0, kind: "none", source: "unsupported" };
  }

  let recoveredImplicit = 0;
  let kind = "none";
  if (foundImplicit) { recoveredImplicit = implicitSum; kind = "implicit"; }
  else if (foundTotal && foundExplicit) { recoveredImplicit = Math.max(0, totalSum - explicitSum); kind = "total_minus_explicit"; }
  else if (foundTotal) { recoveredImplicit = totalSum; kind = "total"; }
  else if (foundExplicit) { recoveredImplicit = explicitSum; kind = "explicit_attached"; }

  return { supported: true, explicitHydrogenCount: explicitSum, implicitHydrogenCount: recoveredImplicit, kind, source: "atom_api" };
}

function parse_json_text_safe(jsonText) {
  if (typeof jsonText !== "string" || !jsonText.trim()) return null;
  try { return JSON.parse(jsonText); } catch (_) { return null; }
}

function normalized_key_map(obj) {
  const out = Object.create(null);
  if (!obj || typeof obj !== "object") return out;
  const keys = Object.keys(obj);
  for (let i = 0; i < keys.length; i++) {
    const raw = keys[i];
    const norm = String(raw).toLowerCase().replace(/[\_\-\s]/g, "");
    if (!(norm in out)) out[norm] = raw;
  }
  return out;
}

function json_record_number(obj, candidateKeys) {
  if (!obj || typeof obj !== "object") return null;
  const keyMap = normalized_key_map(obj);
  for (let i = 0; i < candidateKeys.length; i++) {
    const norm = String(candidateKeys[i]).toLowerCase().replace(/[\_\-\s]/g, "");
    if (!(norm in keyMap)) continue;
    return normalize_count_safe(obj[keyMap[norm]]);
  }
  return null;
}

function looks_like_atom_record(obj) {
  if (!obj || typeof obj !== "object" || Array.isArray(obj)) return false;
  const keyMap = normalized_key_map(obj);
  const keys = Object.keys(keyMap);
  for (let i = 0; i < keys.length; i++) {
    const k = keys[i];
    if (
      k === "atomicnum" || k === "atomicnumber" || k === "symbol" || k === "sym" ||
      k === "element" || k === "el" || k === "atom" || k === "imphs" ||
      k === "numimpliciths" || k === "impliciths" || k === "numexpliciths" ||
      k === "expliciths" || k === "totalhs" || k === "numtotalhs"
    ) return true;
  }
  return false;
}

function get_json_atom_records(root) {
  if (!root || typeof root !== "object") return null;
  const directCandidates = [
    root.atoms,
    root.mol && root.mol.atoms,
    root.molecule && root.molecule.atoms,
    Array.isArray(root.molecules) && root.molecules[0] ? root.molecules[0].atoms : null,
    Array.isArray(root.molecules) && root.molecules[0] && root.molecules[0].mol ? root.molecules[0].mol.atoms : null
  ];
  for (let i = 0; i < directCandidates.length; i++) {
    const arr = directCandidates[i];
    if (Array.isArray(arr) && arr.length) return arr;
  }
  const seen = new Set();
  const queue = [root];
  let guard = 0;
  while (queue.length && guard < 2000) {
    guard += 1;
    const cur = queue.shift();
    if (!cur || typeof cur !== "object") continue;
    if (seen.has(cur)) continue;
    seen.add(cur);
    if (Array.isArray(cur)) {
      if (cur.length && cur.every((it) => looks_like_atom_record(it))) return cur;
      for (let i = 0; i < cur.length; i++) queue.push(cur[i]);
      continue;
    }
    const vals = Object.values(cur);
    for (let i = 0; i < vals.length; i++) queue.push(vals[i]);
  }
  return null;
}

function collect_json_hydrogen_info(mol) {
  const root = parse_json_text_safe(get_json_text(mol));
  const atoms = get_json_atom_records(root);
  if (!Array.isArray(atoms) || !atoms.length) {
    return { supported: false, explicitHydrogenCount: 0, implicitHydrogenCount: 0, kind: "none", source: "unsupported" };
  }

  let explicitSum = 0;
  let implicitSum = 0;
  let totalSum = 0;
  let foundExplicit = false;
  let foundImplicit = false;
  let foundTotal = false;

  for (let i = 0; i < atoms.length; i++) {
    const rec = atoms[i];
    const exp = json_record_number(rec, ["numExplicitHs", "explicitHs", "num_explicit_hs", "nExplicitHs"]);
    const imp = json_record_number(rec, ["impHs", "numImplicitHs", "implicitHs", "num_implicit_hs", "nImplicitHs"]);
    const total = json_record_number(rec, ["totalHs", "numTotalHs", "total_hs"]);
    if (exp != null) { explicitSum += exp; foundExplicit = true; }
    if (imp != null) { implicitSum += imp; foundImplicit = true; }
    if (total != null) { totalSum += total; foundTotal = true; }
  }

  if (!foundExplicit && !foundImplicit && !foundTotal) {
    return { supported: false, explicitHydrogenCount: 0, implicitHydrogenCount: 0, kind: "none", source: "unsupported" };
  }

  let recoveredImplicit = 0;
  let kind = "none";
  if (foundImplicit) { recoveredImplicit = implicitSum; kind = "implicit"; }
  else if (foundTotal && foundExplicit) { recoveredImplicit = Math.max(0, totalSum - explicitSum); kind = "total_minus_explicit"; }
  else if (foundTotal) { recoveredImplicit = totalSum; kind = "total"; }
  else if (foundExplicit) { recoveredImplicit = explicitSum; kind = "explicit_attached"; }

  return { supported: true, explicitHydrogenCount: explicitSum, implicitHydrogenCount: recoveredImplicit, kind, source: "json" };
}

function merge_unique_symbols(baseSymbols, extraSym) {
  const uniq = new Set(Array.isArray(baseSymbols) ? baseSymbols : []);
  if (extraSym) uniq.add(extraSym);
  return Array.from(uniq).sort();
}

function make_identity_stub() {
  return {
    Zs: [],
    syms: [],
    atomCount: 0,
    hydrogenCount: 0,
    uniqueSymbols: [],
    uiUniqueSymbols: [],
    explicitHydrogenCount: 0,
    implicitHydrogenCount: 0,
    hydrogenMode: "unsupported",
    hasHydrogenElement: false,
    explicitHWorked: false,
    source: "fallback",
    hydrogenInfoSource: "unsupported",
    hydrogenInfoKind: "none"
  };
}


function get_mol_identity_info(mol) {
  if (!mol) return null;

  const replacement = replacement_get_cache(mol);
  if (replacement && replacement.identity) return clone_identity_info(replacement.identity);

  const atomCount = get_num_atoms_safe(mol);
  if (!(atomCount > 0)) return make_identity_stub();

  const parsed = get_atoms_bonds_safe(mol);
  const parsedAtoms = Array.isArray(parsed && parsed.atoms) ? parsed.atoms : [];
  const zs = new Array(atomCount);
  const syms = new Array(atomCount);
  let source = "atom_api";

  for (let i = 0; i < atomCount; i++) {
    let z = 0;
    let sym = "";
    const atom = raw_get_atom_with_idx(mol, i);
    if (atom) {
      z = raw_atom_get_atomic_number(atom);
      sym = raw_atom_get_symbol(atom);
      safe_delete(atom);
    }

    if (!(z > 0) || !sym) {
      source = "molblock";
      const rec = parsedAtoms[i] || null;
      if (!(z > 0) && rec && Number.isFinite(rec.Z)) z = rec.Z | 0;
      if (!sym && rec && rec.sym) sym = normalize_symbol(rec.sym);
    }

    if (!(z > 0) && sym) z = symbol_to_Z_safe(sym);
    if (!sym && z > 0 && ELEMENT_SYMBOLS[z]) sym = ELEMENT_SYMBOLS[z];

    zs[i] = z > 0 ? z : 0;
    syms[i] = sym || null;
  }

  const uniq = new Set();
  let explicitHydrogenCount = 0;
  for (let i = 0; i < atomCount; i++) {
    const sym = normalize_symbol(syms[i]);
    const z = zs[i] | 0;
    if (sym) uniq.add(sym);
    if (sym === "H" || z === 1) explicitHydrogenCount += 1;
  }

  const uniqueSymbols = Array.from(uniq).sort();
  let implicitHydrogenCount = 0;
  let implicitHydrogenSupported = false;
  let hydrogenInfoSource = explicitHydrogenCount > 0 ? "explicit_molblock" : "unsupported";
  let hydrogenInfoKind = explicitHydrogenCount > 0 ? "explicit" : "none";

  if (explicitHydrogenCount === 0) {
    const atomApiH = collect_atom_api_hydrogen_info(mol, atomCount);
    if (atomApiH.supported) {
      implicitHydrogenSupported = true;
      implicitHydrogenCount = atomApiH.implicitHydrogenCount;
      hydrogenInfoSource = atomApiH.source;
      hydrogenInfoKind = atomApiH.implicitHydrogenCount > 0 ? atomApiH.kind : "none";
    } else {
      const jsonH = collect_json_hydrogen_info(mol);
      if (jsonH.supported) {
        implicitHydrogenSupported = true;
        implicitHydrogenCount = jsonH.implicitHydrogenCount;
        hydrogenInfoSource = jsonH.source;
        hydrogenInfoKind = jsonH.implicitHydrogenCount > 0 ? jsonH.kind : "none";
      }
    }
  }

  let hydrogenMode = "unsupported";
  let uiUniqueSymbols = uniqueSymbols.slice(0);
  let hydrogenCount = explicitHydrogenCount;
  let hasHydrogenElement = uniqueSymbols.includes("H");

  if (explicitHydrogenCount > 0) {
    hydrogenMode = "explicit";
    hasHydrogenElement = true;
  } else if (implicitHydrogenCount > 0) {
    hydrogenMode = "implicit_only";
    uiUniqueSymbols = merge_unique_symbols(uniqueSymbols, "H");
    hydrogenCount = implicitHydrogenCount;
    hasHydrogenElement = true;
  } else if (implicitHydrogenSupported) {
    hydrogenMode = "none";
  }

  const info = {
    Zs: zs,
    syms: syms,
    atomCount,
    hydrogenCount,
    uniqueSymbols,
    uiUniqueSymbols,
    explicitHydrogenCount,
    implicitHydrogenCount,
    hydrogenMode,
    hasHydrogenElement,
    explicitHWorked: false,
    source,
    hydrogenInfoSource,
    hydrogenInfoKind
  };

  dlog("H policy", {
    hydrogenMode: info.hydrogenMode,
    explicitHydrogenCount: info.explicitHydrogenCount,
    implicitHydrogenCount: info.implicitHydrogenCount,
    uiUniqueSymbols: info.uiUniqueSymbols,
    source: info.hydrogenInfoSource
  });

  return info;
}

function clone_identity_info(info) {
  const src = info && typeof info === "object" ? info : make_identity_stub();
  return {
    Zs: Array.isArray(src.Zs) ? src.Zs.slice(0) : [],
    syms: Array.isArray(src.syms) ? src.syms.slice(0) : [],
    atomCount: Number.isFinite(src.atomCount) ? (src.atomCount | 0) : 0,
    hydrogenCount: Number.isFinite(src.hydrogenCount) ? (src.hydrogenCount | 0) : 0,
    uniqueSymbols: Array.isArray(src.uniqueSymbols) ? src.uniqueSymbols.slice(0) : [],
    uiUniqueSymbols: Array.isArray(src.uiUniqueSymbols) ? src.uiUniqueSymbols.slice(0) : [],
    explicitHydrogenCount: Number.isFinite(src.explicitHydrogenCount) ? (src.explicitHydrogenCount | 0) : 0,
    implicitHydrogenCount: Number.isFinite(src.implicitHydrogenCount) ? (src.implicitHydrogenCount | 0) : 0,
    hydrogenMode: src.hydrogenMode || "unsupported",
    hasHydrogenElement: !!src.hasHydrogenElement,
    explicitHWorked: !!src.explicitHWorked,
    source: src.source || "fallback",
    hydrogenInfoSource: src.hydrogenInfoSource || "unsupported",
    hydrogenInfoKind: src.hydrogenInfoKind || "none"
  };
}

function clone_atom_record(atom) {
  const a = atom && typeof atom === "object" ? atom : {};
  const out = {
    atom_idx: Number.isFinite(a.atom_idx) ? (a.atom_idx | 0) : 0,
    Z: Number.isFinite(a.Z) ? (a.Z | 0) : 0,
    x: Number.isFinite(a.x) ? a.x : Number(a.x) || 0,
    y: Number.isFinite(a.y) ? a.y : Number(a.y) || 0,
    z: Number.isFinite(a.z) ? a.z : Number(a.z) || 0
  };
  if (a.source_atom_id != null) out.source_atom_id = Number.isFinite(a.source_atom_id) ? (a.source_atom_id | 0) : parse_int_safe(a.source_atom_id, 0);
  if (a.sym) out.sym = normalize_symbol(a.sym);
  return out;
}

function clone_coord_record(coord) {
  const c = coord && typeof coord === "object" ? coord : {};
  return {
    x: Number.isFinite(c.x) ? c.x : Number(c.x) || 0,
    y: Number.isFinite(c.y) ? c.y : Number(c.y) || 0,
    z: Number.isFinite(c.z) ? c.z : Number(c.z) || 0
  };
}

function clone_bond_record(bond) {
  if (!Array.isArray(bond)) return null;
  return [parse_int_safe(bond[0], 0), parse_int_safe(bond[1], 0), parse_bond_order_token(bond[2])];
}

function push_unique_warning(warnings, msg) {
  const text = String(msg == null ? "" : msg).trim();
  if (!text) return;
  if (!Array.isArray(warnings)) return;
  if (warnings.indexOf(text) >= 0) return;
  warnings.push(text);
}

function make_warning_groups() {
  return { input: [], backend: [], geometry: [] };
}

function clone_warning_groups(groups) {
  const src = groups && typeof groups === "object" ? groups : {};
  return {
    input: Array.isArray(src.input) ? src.input.slice(0) : [],
    backend: Array.isArray(src.backend) ? src.backend.slice(0) : [],
    geometry: Array.isArray(src.geometry) ? src.geometry.slice(0) : []
  };
}

function push_group_warning(groups, warnings, groupName, msg) {
  push_unique_warning(warnings, msg);
  if (!groups || typeof groups !== "object") return;
  if (!Array.isArray(groups[groupName])) groups[groupName] = [];
  push_unique_warning(groups[groupName], msg);
}

function count_explicit_h_atoms_in_atoms(atoms) {
  let n = 0;
  const arr = Array.isArray(atoms) ? atoms : [];
  for (let i = 0; i < arr.length; i++) {
    const sym = normalize_symbol(arr[i] && arr[i].sym);
    const z = arr[i] && Number.isFinite(arr[i].Z) ? (arr[i].Z | 0) : parse_int_safe(arr[i] && arr[i].Z, 0);
    if (sym === "H" || z === 1) n += 1;
  }
  return n;
}

function count_symbol_z_mismatches(identity) {
  const info = identity && typeof identity === "object" ? identity : {};
  const zs = Array.isArray(info.Zs) ? info.Zs : [];
  const syms = Array.isArray(info.syms) ? info.syms : [];
  const n = Math.max(zs.length, syms.length);
  let missingZForSymbol = 0;
  let missingSymForZ = 0;
  for (let i = 0; i < n; i++) {
    const z = Number.isFinite(zs[i]) ? (zs[i] | 0) : parse_int_safe(zs[i], 0);
    const sym = normalize_symbol(syms[i]);
    if (sym && !(z > 0)) missingZForSymbol += 1;
    if (!sym && z > 0) missingSymForZ += 1;
  }
  return { missingZForSymbol, missingSymForZ };
}

function count_unknown_zero_Z_atoms(atoms) {
  let n = 0;
  const arr = Array.isArray(atoms) ? atoms : [];
  for (let i = 0; i < arr.length; i++) {
    const z = arr[i] && Number.isFinite(arr[i].Z) ? (arr[i].Z | 0) : parse_int_safe(arr[i] && arr[i].Z, 0);
    if (!(z > 0)) n += 1;
  }
  return n;
}

function make_chemistry_assessment_stub() {
  return {
    topologyTrusted: false,
    identityTrusted: false,
    hydrogenTrusted: false,
    geometryTrusted: false,
    renderTrusted: false,
    strictOk: false,
    trustLevel: "rejected",
    rejectReasons: []
  };
}

function clone_chemistry_assessment(chemistry) {
  const src = chemistry && typeof chemistry === "object" ? chemistry : make_chemistry_assessment_stub();
  return {
    topologyTrusted: !!src.topologyTrusted,
    identityTrusted: !!src.identityTrusted,
    hydrogenTrusted: !!src.hydrogenTrusted,
    geometryTrusted: !!src.geometryTrusted,
    renderTrusted: !!src.renderTrusted,
    strictOk: !!src.strictOk,
    trustLevel: src.trustLevel || "rejected",
    rejectReasons: Array.isArray(src.rejectReasons) ? src.rejectReasons.slice(0) : []
  };
}

function has_warning(list, code) {
  return Array.isArray(list) && list.indexOf(code) >= 0;
}

function push_reject_reason(out, code) {
  if (!out || !Array.isArray(out.rejectReasons)) return;
  push_unique_warning(out.rejectReasons, code);
}

function resolve_backend_from_input_for_trust(input) {
  if (!input || typeof input !== "object") return null;
  if (input.backend && input.backend.mol) return input.backend;
  if (input.__tem_smiles_backend && input.__tem_smiles_backend.mol) return input.__tem_smiles_backend;
  if (input.mol && input.identity && input.stages) return input;
  return null;
}

function resolve_geometry_from_input_for_trust(input) {
  if (!input || typeof input !== "object") return null;
  if (input.geometry && (input.geometry.mol || typeof input.geometry.geometryMode === "string")) return input.geometry;
  if (input.__tem_smiles_geometry && (input.__tem_smiles_geometry.mol || typeof input.__tem_smiles_geometry.geometryMode === "string")) return input.__tem_smiles_geometry;
  if (typeof input.geometryMode === "string" || typeof input.coordsAvailable === "boolean") return input;
  return null;
}


function assess_smiles_backend_trust(result) {
  const out = make_chemistry_assessment_stub();
  const backend = result && typeof result === "object" ? result : null;
  if (!backend) {
    push_reject_reason(out, "backend_invalid");
    return out;
  }

  const warnings = Array.isArray(backend.warnings) ? backend.warnings : [];
  const atoms = Array.isArray(backend.atoms) ? backend.atoms : [];
  const identity = backend.identity && typeof backend.identity === "object" ? backend.identity : make_identity_stub();
  const unknownZeroZCount = count_unknown_zero_Z_atoms(atoms);

  const topologyTrusted = !!(
    backend.ok &&
    !has_warning(warnings, "mol_create_failed") &&
    !has_warning(warnings, "mol_has_no_atoms") &&
    !has_warning(warnings, "atom_count_mismatch") &&
    !has_warning(warnings, "bond_count_mismatch") &&
    !has_warning(warnings, "replacement_heavy_atom_mismatch_after_bridge") &&
    !has_warning(warnings, "replacement_heavy_atom_mismatch_after_materialization") &&
    !has_warning(warnings, "replacement_heavy_atom_count_mismatch") &&
    !has_warning(warnings, "replacement_heavy_atom_symbol_mismatch") &&
    !has_warning(warnings, "replacement_heavy_bond_graph_mismatch") &&
    !has_warning(warnings, "replacement_heavy_local_signature_mismatch")
  );

  const identityTrusted = !!(
    topologyTrusted &&
    !has_warning(warnings, "symbol_without_Z") &&
    !has_warning(warnings, "Z_without_symbol") &&
    !(unknownZeroZCount > 0)
  );

  const hydrogenTrusted = !!(
    identityTrusted &&
    (identity.hydrogenMode === "explicit" || identity.hydrogenMode === "none") &&
    !has_warning(warnings, "explicit_h_reported_but_not_materialized") &&
    !has_warning(warnings, "explicit_h_mode_without_h_atoms") &&
    !has_warning(warnings, "implicit_only_mode_but_h_atoms_present") &&
    !has_warning(warnings, "hydrogen_explicit_implicit_mixed") &&
    !has_warning(warnings, "replacement_explicit_h_materialization_incomplete")
  );

  out.topologyTrusted = topologyTrusted;
  out.identityTrusted = identityTrusted;
  out.hydrogenTrusted = hydrogenTrusted;
  out.geometryTrusted = false;
  out.renderTrusted = false;
  out.strictOk = false;

  if (has_warning(warnings, "mol_create_failed")) push_reject_reason(out, "mol_create_failed");
  if (has_warning(warnings, "mol_has_no_atoms")) push_reject_reason(out, "mol_has_no_atoms");
  if (has_warning(warnings, "atom_count_mismatch")) push_reject_reason(out, "atom_count_mismatch");
  if (has_warning(warnings, "bond_count_mismatch")) push_reject_reason(out, "bond_count_mismatch");
  if (has_warning(warnings, "symbol_without_Z")) push_reject_reason(out, "identity_symbol_without_Z");
  if (has_warning(warnings, "Z_without_symbol")) push_reject_reason(out, "identity_Z_without_symbol");
  if (unknownZeroZCount > 0) push_reject_reason(out, "identity_unknown_zero_Z");
  if (!(identity.hydrogenMode === "explicit" || identity.hydrogenMode === "none")) push_reject_reason(out, "hydrogen_not_explicit");
  if (has_warning(warnings, "explicit_h_reported_but_not_materialized")) push_reject_reason(out, "hydrogen_explicit_not_materialized");
  if (has_warning(warnings, "explicit_h_mode_without_h_atoms")) push_reject_reason(out, "hydrogen_explicit_mode_without_h_atoms");
  if (has_warning(warnings, "implicit_only_mode_but_h_atoms_present")) push_reject_reason(out, "hydrogen_mode_inconsistent");
  if (has_warning(warnings, "hydrogen_explicit_implicit_mixed")) push_reject_reason(out, "hydrogen_explicit_implicit_mixed");
  if (has_warning(warnings, "replacement_explicit_h_materialization_unavailable")) push_reject_reason(out, "replacement_explicit_h_materialization_unavailable");
  if (has_warning(warnings, "replacement_explicit_h_bridge_failed")) push_reject_reason(out, "replacement_explicit_h_bridge_failed");
  if (has_warning(warnings, "replacement_explicit_h_roundtrip_failed")) push_reject_reason(out, "replacement_explicit_h_roundtrip_failed");
  if (has_warning(warnings, "replacement_heavy_atom_mismatch_after_bridge")) push_reject_reason(out, "replacement_heavy_atom_mismatch_after_bridge");
  if (has_warning(warnings, "replacement_heavy_atom_mismatch_after_materialization")) push_reject_reason(out, "replacement_heavy_atom_mismatch_after_materialization");
  if (has_warning(warnings, "replacement_geometry_regressed_after_materialization")) push_reject_reason(out, "replacement_geometry_regressed_after_materialization");

  if (!topologyTrusted) out.trustLevel = "rejected";
  else if (!identityTrusted && !hydrogenTrusted) out.trustLevel = "topology_only";
  else if (!identityTrusted) out.trustLevel = "identity_degraded";
  else if (!hydrogenTrusted) out.trustLevel = "hydrogen_degraded";
  else out.trustLevel = "geometry_degraded";

  return out;
}



function make_empty_smiles_backend_result(inputSmiles, normalizedSmiles, token) {
  return {
    ok: false,
    smiles: token,
    inputSmiles,
    normalizedSmiles,
    mol: null,
    atoms: [],
    bonds: [],
    coords: [],
    geometryOrigin: "failed",
    authoritativeTopologySource: "fallback",
    authoritativeHydrogenSource: "fallback",
    authoritativeGeometrySource: "fallback",
    authoritativeMolKind: "fallback",
    identity: make_identity_stub(),
    stages: {
      molCreated: false,
      addHsAttempted: false,
      addHsChanged: false,
      coordsAvailable: false,
      atomsAvailable: false,
      bondsAvailable: false
    },
    diagnostics: {
      atomCountFromMol: 0,
      atomCountFromAtoms: 0,
      bondCountFromMol: 0,
      bondCountFromBonds: 0,
      coordsCount: 0,
      molblockSummary: { atomCount: 0, hydrogenCount: 0, uniqueSymbols: [], format: "none" },
      source: "fallback",
      authoritativeTopologySource: "fallback",
      authoritativeHydrogenSource: "fallback",
      authoritativeGeometrySource: "fallback",
      authoritativeMolKind: "fallback",
      symbolZMismatchCount: 0,
      missingZForSymbolCount: 0,
      missingSymForZCount: 0,
      unknownZeroZCount: 0
    },
    warnings: [],
    warningGroups: make_warning_groups(),
    chemistry: make_chemistry_assessment_stub()
  };
}




function build_smiles_backend_result(smiles) {
  const inputSmiles = String(smiles == null ? "" : smiles);
  const normalizedSmiles = inputSmiles.trim();
  const token = normalizedSmiles || "O";
  const result = make_empty_smiles_backend_result(inputSmiles, normalizedSmiles, token);

  if (!normalizedSmiles) {
    push_group_warning(result.warningGroups, result.warnings, "input", "empty_input_fallback_to_O");
  }

  let replacement = replacement_extract_smiles_result(token);

  let mol = raw_get_mol(token);
  if (!mol && replacement && replacement.ok) {
    mol = make_adapter_mol_wrapper(token);
    replacement.usedWrapperMol = true;
    push_group_warning(result.warningGroups, result.warnings, "backend", "rdkit_mol_unavailable_using_wrapper");
  }

  if (!mol) {
    push_group_warning(result.warningGroups, result.warnings, "backend", "mol_create_failed");
    result.chemistry = assess_smiles_backend_trust(result);
    dlog("SMILES backend result", { ok: result.ok, warnings: result.warnings });
    return result;
  }

  result.stages.molCreated = true;
  result.mol = mol;

  if (mol) {
    const replacementNeedsExplicitHydrogen =
      replacement &&
      replacement.ok &&
      replacement.identity &&
      replacement.identity.hydrogenMode !== "explicit" &&
      replacement.identity.hydrogenMode !== "none";
    if ((!replacement || !replacement.ok || replacementNeedsExplicitHydrogen)) {
      const directFallback = build_rdkit_explicit_h_fallback_payload(mol);
      if (directFallback && directFallback.payload && directFallback.payload.ok) {
        if (directFallback.mol && directFallback.mol !== mol) {
          safe_delete(mol);
          mol = directFallback.mol;
          result.mol = mol;
        }
        replacement = directFallback.payload;
      }
    }
  }

  let hsRes = null;
  if (replacement && replacement.ok) {
    replacement_attach_to_mol(mol, replacement);
    hsRes = {
      mol,
      changed: !!replacement.expandWorked,
      apiAvailable: false,
      returnedNewMol: false,
      atomCountBefore: replacement.baseAtomCount | 0,
      atomCountAfter: Array.isArray(replacement.atoms) ? replacement.atoms.length : 0,
      molblockHydrogenCountAfter: replacement.explicitHydrogenCount | 0,
      source: replacement.authoritativeHydrogenSource || replacement.backendName || REPLACEMENT_BACKEND_NAME
    };
    result.stages.addHsAttempted = !!replacement.expandWorked;
    result.stages.addHsChanged = !!replacement.expandWorked;
  } else {
    hsRes = add_hs_safe(mol);
    result.stages.addHsAttempted = !!(hsRes && hsRes.apiAvailable);
    result.stages.addHsChanged = !!(hsRes && hsRes.changed);

    const outMol = hsRes && hsRes.mol ? hsRes.mol : mol;
    if (outMol && outMol !== mol) {
      safe_delete(mol);
      mol = outMol;
      result.mol = mol;
    }
  }

  const identity = replacement && replacement.identity
    ? clone_identity_info(replacement.identity)
    : clone_identity_info(get_mol_identity_info(mol));

  identity.explicitHWorked = replacement && replacement.identity
    ? !!replacement.identity.explicitHWorked
    : !!(
        hsRes &&
        (hsRes.molblockHydrogenCountAfter > 0 || hsRes.atomCountAfter > hsRes.atomCountBefore)
      );
  result.identity = identity;

  const atomsRaw = replacement && replacement.ok
    ? {
        atoms: Array.isArray(replacement.atoms) ? replacement.atoms.map(clone_atom_record) : [],
        bonds: Array.isArray(replacement.bonds) ? replacement.bonds.map(clone_bond_record).filter(Boolean) : []
      }
    : get_atoms_bonds_safe(mol);

  const atoms = Array.isArray(atomsRaw && atomsRaw.atoms) ? atomsRaw.atoms.map(clone_atom_record) : [];
  const bonds = Array.isArray(atomsRaw && atomsRaw.bonds) ? atomsRaw.bonds.map(clone_bond_record).filter(Boolean) : [];

  const rdkitCoords = get_molblock_coordinates_only(mol).map(clone_coord_record);
  const rdkitCoordsUsable = rdkitCoords.length === atoms.length && rdkitCoords.length > 0 && !are_atoms_degenerate(rdkitCoords);
  const replacementCoords = replacement && replacement.ok && Array.isArray(replacement.coords)
    ? replacement.coords.map(clone_coord_record)
    : [];
  const replacementCoordsUsable = replacementCoords.length === atoms.length && replacementCoords.length > 0 && !are_atoms_degenerate(replacementCoords);

  dlog("SMILES backend coord candidates", {
    atomCount: atoms.length,
    rdkitCoordCount: rdkitCoords.length,
    rdkitCoordsUsable,
    replacementCoordCount: replacementCoords.length,
    replacementCoordsUsable,
    preferReplacementCoords: replacementCoordsUsable && (!rdkitCoordsUsable || replacementCoords.length > rdkitCoords.length)
  });

  let coords = [];
  let geometryOrigin = "rdkit";
  if (replacementCoordsUsable && (!rdkitCoordsUsable || replacementCoords.length > rdkitCoords.length)) {
    coords = replacementCoords;
    geometryOrigin = replacement.geometryOrigin || "replacement_2d";
  } else if (rdkitCoordsUsable) {
    coords = rdkitCoords;
    geometryOrigin = "rdkit";
  } else if (replacementCoordsUsable) {
    coords = replacementCoords;
    geometryOrigin = replacement.geometryOrigin || "replacement_2d";
  } else if (rdkitCoords.length === atoms.length && rdkitCoords.length > 0) {
    coords = rdkitCoords;
    geometryOrigin = "rdkit";
  } else if (replacementCoords.length === atoms.length && replacementCoords.length > 0) {
    coords = replacementCoords;
    geometryOrigin = replacement.geometryOrigin || "replacement_2d";
  } else {
    coords = get_coordinates_safe(mol).map(clone_coord_record);
    geometryOrigin = replacement && replacement.ok ? (replacement.geometryOrigin || "replacement_none") : "rdkit";
  }

  result.atoms = atoms;
  result.bonds = bonds;
  result.coords = coords;
  result.geometryOrigin = geometryOrigin;
  result.authoritativeTopologySource = replacement && replacement.ok
    ? (replacement.authoritativeTopologySource || "openchemlib")
    : "rdkit";
  result.authoritativeHydrogenSource = replacement && replacement.ok
    ? (replacement.authoritativeHydrogenSource || (identity.hydrogenMode === "explicit" ? "openchemlib" : "fallback"))
    : (identity.hydrogenMode === "explicit" ? "rdkit" : "fallback");
  result.authoritativeGeometrySource = rdkitCoordsUsable
    ? "rdkit"
    : replacement && replacement.ok
      ? (replacement.authoritativeGeometrySource || "fallback")
      : "rdkit";
  result.authoritativeMolKind = replacement_is_wrapper_mol(mol)
    ? "wrapper_only"
    : replacement && replacement.ok
      ? (replacement.authoritativeMolKind || "rdkit_token_mol")
      : "rdkit_token_mol";

  result.stages.atomsAvailable = atoms.length > 0;
  result.stages.bondsAvailable = bonds.length > 0;
  result.stages.coordsAvailable = coords.length > 0;

  const atomCountFromMol = replacement && replacement.ok ? atoms.length : get_num_atoms_safe(mol);
  const bondCountFromMol = replacement && replacement.ok ? bonds.length : get_num_bonds_safe(mol);
  const molblockSummary = summarize_molblock(mol);
  const mismatch = count_symbol_z_mismatches(identity);
  const explicitHydrogenCountFromAtoms = count_explicit_h_atoms_in_atoms(atoms);
  const unknownZeroZCount = count_unknown_zero_Z_atoms(atoms);

  result.diagnostics = {
    atomCountFromMol,
    atomCountFromAtoms: atoms.length,
    bondCountFromMol,
    bondCountFromBonds: bonds.length,
    coordsCount: coords.length,
    molblockSummary,
    source: replacement && replacement.ok ? (replacement.source || "replacement_backend") : (identity.source || "fallback"),
    authoritativeTopologySource: result.authoritativeTopologySource,
    authoritativeHydrogenSource: result.authoritativeHydrogenSource,
    authoritativeGeometrySource: result.authoritativeGeometrySource,
    authoritativeMolKind: result.authoritativeMolKind,
    symbolZMismatchCount: mismatch.missingZForSymbol + mismatch.missingSymForZ,
    missingZForSymbolCount: mismatch.missingZForSymbol,
    missingSymForZCount: mismatch.missingSymForZ,
    unknownZeroZCount,
    replacementBackend: replacement && replacement.ok ? {
      name: replacement.backendName || REPLACEMENT_BACKEND_NAME,
      source: replacement_get_status().source || "preloaded",
      usedWrapperMol: !!replacement.usedWrapperMol,
      baseAtomCount: replacement.baseAtomCount | 0,
      baseBondCount: replacement.baseBondCount | 0,
      hydrogenAddedCount: replacement.hydrogenAddedCount | 0,
      authoritativeTopologySource: replacement.authoritativeTopologySource || "openchemlib",
      authoritativeHydrogenSource: replacement.authoritativeHydrogenSource || "fallback",
      authoritativeGeometrySource: replacement.authoritativeGeometrySource || "fallback",
      authoritativeMolKind: replacement.authoritativeMolKind || "fallback",
      bridgeRepresentation: replacement.diagnostics && replacement.diagnostics.bridgeRepresentation ? replacement.diagnostics.bridgeRepresentation : ""
    } : null,
    rdkitRawAtomCount: raw_get_num_atoms(mol),
    rdkitRawBondCount: raw_get_num_bonds(mol)
  };

  if (!(atomCountFromMol > 0)) push_group_warning(result.warningGroups, result.warnings, "backend", "mol_has_no_atoms");
  if (!(replacement && replacement.ok) && atomCountFromMol !== atoms.length) push_group_warning(result.warningGroups, result.warnings, "backend", "atom_count_mismatch");
  if (!(replacement && replacement.ok) && bondCountFromMol !== bonds.length) push_group_warning(result.warningGroups, result.warnings, "backend", "bond_count_mismatch");
  if (atomCountFromMol > 0 && coords.length === 0) push_group_warning(result.warningGroups, result.warnings, "backend", "coords_missing");
  if (coords.length > 0 && atoms.length > 0 && coords.length !== atoms.length) push_group_warning(result.warningGroups, result.warnings, "backend", "coords_atom_count_mismatch");
  if (mismatch.missingZForSymbol > 0) push_group_warning(result.warningGroups, result.warnings, "backend", "symbol_without_Z");
  if (mismatch.missingSymForZ > 0) push_group_warning(result.warningGroups, result.warnings, "backend", "Z_without_symbol");
  if (identity.explicitHWorked && explicitHydrogenCountFromAtoms === 0 && identity.hydrogenCount > 0) push_group_warning(result.warningGroups, result.warnings, "backend", "explicit_h_reported_but_not_materialized");
  if (identity.hydrogenMode === "explicit" && explicitHydrogenCountFromAtoms === 0 && identity.hydrogenCount > 0) push_group_warning(result.warningGroups, result.warnings, "backend", "explicit_h_mode_without_h_atoms");
  if (identity.hydrogenMode === "implicit_only" && explicitHydrogenCountFromAtoms > 0) push_group_warning(result.warningGroups, result.warnings, "backend", "implicit_only_mode_but_h_atoms_present");
  if (identity.hydrogenMode === "mixed") push_group_warning(result.warningGroups, result.warnings, "backend", "hydrogen_explicit_implicit_mixed");
  if (!(replacement && replacement.ok) && !result.stages.addHsAttempted) push_group_warning(result.warningGroups, result.warnings, "backend", "addHs_not_available");
  if (replacement && replacement.ok && Array.isArray(replacement.warnings)) {
    for (let i = 0; i < replacement.warnings.length; i++) {
      push_group_warning(result.warningGroups, result.warnings, "backend", replacement.warnings[i]);
    }
  }

  result.ok = !!(result.stages.molCreated && identity.atomCount > 0);
  result.chemistry = assess_smiles_backend_trust(result);

  dlog("SMILES backend result", {
    ok: result.ok,
    atomCountFromMol: result.diagnostics.atomCountFromMol,
    atomCountFromAtoms: result.diagnostics.atomCountFromAtoms,
    bondCountFromMol: result.diagnostics.bondCountFromMol,
    bondCountFromBonds: result.diagnostics.bondCountFromBonds,
    coordsAvailable: result.stages.coordsAvailable,
    hydrogenMode: result.identity.hydrogenMode,
    uniqueSymbols: result.identity.uniqueSymbols,
    uiUniqueSymbols: result.identity.uiUniqueSymbols,
    geometryOrigin: result.geometryOrigin,
    authoritativeTopologySource: result.authoritativeTopologySource,
    authoritativeHydrogenSource: result.authoritativeHydrogenSource,
    authoritativeGeometrySource: result.authoritativeGeometrySource,
    authoritativeMolKind: result.authoritativeMolKind,
    replacementBackend: result.diagnostics.replacementBackend,
    warnings: result.warnings
  });

  return result;
}


function validate_smiles_backend(smiles) {
  const result = build_smiles_backend_result(smiles);
  try {
    return {
      ok: !!(result && result.ok),
      reason:
        result && result.ok
          ? "ok"
          : result && result.chemistry && Array.isArray(result.chemistry.rejectReasons) && result.chemistry.rejectReasons.length
            ? result.chemistry.rejectReasons[0]
            : result && Array.isArray(result.warnings) && result.warnings.length
              ? result.warnings[0]
              : "backend_invalid",
      atomCount: result && result.identity ? (result.identity.atomCount | 0) : 0,
      warnings: result && Array.isArray(result.warnings) ? result.warnings.slice(0) : [],
      rejectReasons: result && result.chemistry ? result.chemistry.rejectReasons.slice(0) : [],
      trustLevel: result && result.chemistry ? result.chemistry.trustLevel : "rejected"
    };
  } finally {
    if (result && result.mol) safe_delete(result.mol);
  }
}

function has_usable_smiles_geometry(result) {
  const r = result && typeof result === "object" ? result : null;
  if (!r) return false;
  if (typeof r.coordsAvailable === "boolean") {
    const atoms0 = Array.isArray(r.atoms) ? r.atoms : [];
    const coords0 = Array.isArray(r.coords) ? r.coords : [];
    if (!r.coordsAvailable) return false;
    if (!(atoms0.length > 0) || !(coords0.length > 0)) return false;
    return atoms0.length === coords0.length;
  }
  if (!r.ok) return false;
  if (!r.stages || !r.stages.coordsAvailable || !r.stages.atomsAvailable) return false;
  const atoms = Array.isArray(r.atoms) ? r.atoms : [];
  const coords = Array.isArray(r.coords) ? r.coords : [];
  if (!(atoms.length > 0) || !(coords.length > 0)) return false;
  return atoms.length === coords.length;
}

function atoms_from_identity_and_coords(identity, coords) {
  const info = identity && typeof identity === "object" ? identity : make_identity_stub();
  const out = [];
  const arr = Array.isArray(coords) ? coords : [];
  for (let i = 0; i < arr.length; i++) {
    const c = clone_coord_record(arr[i]);
    const z = Array.isArray(info.Zs) && Number.isFinite(info.Zs[i]) ? (info.Zs[i] | 0) : 0;
    const sym = Array.isArray(info.syms) ? normalize_symbol(info.syms[i]) : "";
    const atom = { atom_idx: i, Z: z, x: c.x, y: c.y, z: c.z };
    if (sym) atom.sym = sym;
    out.push(atom);
  }
  return out;
}

function atoms_from_identity_circle(identity, radius = 3.0) {
  const info = identity && typeof identity === "object" ? identity : make_identity_stub();
  const n = Math.max(0, Number.isFinite(info.atomCount) ? (info.atomCount | 0) : 0);
  const out = [];
  if (!(n > 0)) return out;
  for (let i = 0; i < n; i++) {
    const t = (2 * Math.PI * i) / n;
    const z = Array.isArray(info.Zs) && Number.isFinite(info.Zs[i]) ? (info.Zs[i] | 0) : 0;
    const sym = Array.isArray(info.syms) ? normalize_symbol(info.syms[i]) : "";
    const atom = { atom_idx: i, Z: z, x: Math.cos(t) * radius, y: Math.sin(t) * radius, z: 0 };
    if (sym) atom.sym = sym;
    out.push(atom);
  }
  return out;
}

function are_atoms_degenerate(atoms, eps = 1e-6) {
  const arr = Array.isArray(atoms) ? atoms : [];
  if (arr.length < 2) return arr.length > 0;
  let minx = Number(arr[0] && arr[0].x) || 0;
  let maxx = minx;
  let miny = Number(arr[0] && arr[0].y) || 0;
  let maxy = miny;
  for (let i = 1; i < arr.length; i++) {
    const a = arr[i] || {};
    const x = Number(a.x);
    const y = Number(a.y);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return true;
    if (x < minx) minx = x;
    if (x > maxx) maxx = x;
    if (y < miny) miny = y;
    if (y > maxy) maxy = y;
  }
  return maxx - minx < eps && maxy - miny < eps;
}

function make_empty_smiles_geometry_result(mol) {
  return {
    ok: false,
    mol: mol || null,
    backend: null,
    atoms: [],
    bonds: [],
    coords: [],
    geometryMode: "failed",
    coordsAvailable: false,
    degenerateCoords: false,
    usedFallbackCoords: false,
    warnings: [],
    warningGroups: make_warning_groups(),
    chemistry: make_chemistry_assessment_stub()
  };
}

function refresh_geometry_snapshot(mol, identity, prevBonds) {
  const extracted = get_atoms_bonds_safe(mol);
  let atoms = Array.isArray(extracted && extracted.atoms) ? extracted.atoms.map(clone_atom_record) : [];
  let bonds = Array.isArray(extracted && extracted.bonds) ? extracted.bonds.map(clone_bond_record).filter(Boolean) : [];
  const coords = get_coordinates_safe(mol).map(clone_coord_record);
  if ((!atoms || !atoms.length) && coords.length) atoms = atoms_from_identity_and_coords(identity, coords);
  if ((!bonds || !bonds.length) && Array.isArray(prevBonds) && prevBonds.length) {
    bonds = prevBonds.map(clone_bond_record).filter(Boolean);
  }
  return { atoms, bonds, coords };
}

function build_smiles_geometry_result(input) {
  let backend = null;
  if (typeof input === "string") backend = build_smiles_backend_result(input);
  else if (input && typeof input === "object" && input.mol && input.identity) backend = input;
  else if (input && typeof input === "object" && input.__tem_smiles_backend) backend = input.__tem_smiles_backend;
  else if (input && typeof input === "object") {
    backend = {
      ok: true,
      mol: input,
      atoms: [],
      bonds: [],
      coords: [],
      identity: clone_identity_info(get_mol_identity_info(input)),
      warnings: [],
      warningGroups: make_warning_groups(),
      chemistry: make_chemistry_assessment_stub(),
      stages: { molCreated: true, addHsAttempted: false, addHsChanged: false, coordsAvailable: false, atomsAvailable: false, bondsAvailable: false },
      diagnostics: { atomCountFromMol: get_num_atoms_safe(input), atomCountFromAtoms: 0, bondCountFromMol: get_num_bonds_safe(input), bondCountFromBonds: 0, coordsCount: 0, molblockSummary: summarize_molblock(input), source: "fallback", symbolZMismatchCount: 0, missingZForSymbolCount: 0, missingSymForZCount: 0, unknownZeroZCount: 0 }
    };
  }

  const mol = backend && backend.mol ? backend.mol : null;
  const identity = clone_identity_info(backend && backend.identity);
  const out = make_empty_smiles_geometry_result(mol);
  out.backend = backend || null;

  if (!backend || !backend.ok || !mol) {
    push_group_warning(out.warningGroups, out.warnings, "geometry", "geometry_backend_invalid");
    out.chemistry = assess_smiles_trust({ backend, geometry: out });
    return out;
  }

  let atoms = Array.isArray(backend.atoms) ? backend.atoms.map(clone_atom_record) : [];
  let bonds = Array.isArray(backend.bonds) ? backend.bonds.map(clone_bond_record).filter(Boolean) : [];
  let coords = Array.isArray(backend.coords) ? backend.coords.map(clone_coord_record) : [];

  const replacementCache = mol ? replacement_get_cache(mol) : null;
  const replacementCacheCoords =
    replacementCache && Array.isArray(replacementCache.coords)
      ? replacementCache.coords.map(clone_coord_record)
      : [];
  const replacementCacheUsable =
    replacementCacheCoords.length === atoms.length &&
    replacementCacheCoords.length > 0 &&
    !are_atoms_degenerate(replacementCacheCoords);

  if (replacementCacheUsable && coords.length !== replacementCacheCoords.length) {
    coords = replacementCacheCoords;
    dlog("SMILES geometry adopted replacement cache", {
      atomCount: atoms.length,
      coordCount: coords.length,
      geometryOrigin: replacementCache.geometryOrigin || "replacement_2d"
    });
  }

  if ((!atoms || !atoms.length) && coords.length) atoms = atoms_from_identity_and_coords(identity, coords);
  if ((!bonds || !bonds.length) && backend && Array.isArray(backend.bonds) && backend.bonds.length) {
    bonds = backend.bonds.map(clone_bond_record).filter(Boolean);
  }

  let geometryMode = "failed";
  let degenerateCoords = false;
  let usedFallbackCoords = false;

  function adopt(snapshot, mode) {
    atoms = Array.isArray(snapshot.atoms) ? snapshot.atoms : [];
    bonds = Array.isArray(snapshot.bonds) ? snapshot.bonds : [];
    coords = Array.isArray(snapshot.coords) ? snapshot.coords : [];
    geometryMode = mode;
    degenerateCoords = are_atoms_degenerate(coords);
  }

  const nativeUsable = atoms.length > 0 && coords.length > 0 && atoms.length === coords.length;
  const backendGeometryOrigin = backend && typeof backend.geometryOrigin === "string" ? backend.geometryOrigin : "rdkit";
  dlog("SMILES geometry initial candidate", {
    atomCount: atoms.length,
    coordCount: coords.length,
    nativeUsable,
    backendGeometryOrigin,
    authoritativeGeometrySource: backend && backend.authoritativeGeometrySource ? backend.authoritativeGeometrySource : "fallback"
  });
  if (nativeUsable) {
    let adoptedMode = "native";
    if (backendGeometryOrigin === "replacement_2d") adoptedMode = "replacement_2d";
    else if (backendGeometryOrigin === "2d_generated") adoptedMode = "2d_generated";
    else if (
      backendGeometryOrigin === "replacement_3d_conformer" ||
      backendGeometryOrigin === "replacement_3d_mmff94" ||
      backendGeometryOrigin === "replacement_3d_existing" ||
      backendGeometryOrigin === "embed3d" ||
      backendGeometryOrigin === "rdkit"
    ) adoptedMode = backendGeometryOrigin === "rdkit" ? "native" : "embed3d";
    adopt({ atoms, bonds, coords }, adoptedMode);
  }

  if (geometryMode === "failed" || degenerateCoords) {
    if (degenerateCoords) push_group_warning(out.warningGroups, out.warnings, "geometry", "degenerate_coords");
    const embedded = raw_embed(mol);
    if (embedded) raw_uff_optimize(mol);
    if (embedded) {
      const snapshot = refresh_geometry_snapshot(mol, identity, bonds);
      if (snapshot.atoms.length > 0 && snapshot.coords.length > 0 && snapshot.atoms.length === snapshot.coords.length && !are_atoms_degenerate(snapshot.coords)) {
        adopt(snapshot, "embed3d");
      }
    }
  }

  if (geometryMode === "failed" || degenerateCoords) {
    const did2d = raw_compute2d(mol) || raw_generate_aligned_coords(mol);
    if (did2d) {
      const snapshot = refresh_geometry_snapshot(mol, identity, bonds);
      if (snapshot.atoms.length > 0 && snapshot.coords.length > 0 && snapshot.atoms.length === snapshot.coords.length && !are_atoms_degenerate(snapshot.coords)) {
        adopt(snapshot, "2d_generated");
      } else if (snapshot.atoms.length > 0 && snapshot.coords.length > 0 && snapshot.atoms.length === snapshot.coords.length) {
        atoms = snapshot.atoms;
        bonds = snapshot.bonds;
        coords = snapshot.coords;
        degenerateCoords = true;
        push_group_warning(out.warningGroups, out.warnings, "geometry", "degenerate_coords");
      }
    } else {
      push_group_warning(out.warningGroups, out.warnings, "geometry", "geometry_2d_generation_failed");
    }
  }

  if (geometryMode === "failed") {
    if (!(atoms && atoms.length)) atoms = atoms_from_identity_circle(identity, 3.0);
    if (atoms.length > 0) {
      if (!(bonds && bonds.length) && Array.isArray(backend.bonds) && backend.bonds.length) {
        bonds = backend.bonds.map(clone_bond_record).filter(Boolean);
      }
      coords = atoms.map((a) => clone_coord_record(a));
      geometryMode = "fallback_circle";
      usedFallbackCoords = true;
      degenerateCoords = false;
      push_group_warning(out.warningGroups, out.warnings, "geometry", "fallback_2d_layout_applied");
    }
  }

  out.atoms = Array.isArray(atoms) ? atoms : [];
  out.bonds = Array.isArray(bonds) ? bonds : [];
  out.coords = Array.isArray(coords) ? coords : [];
  out.geometryMode = geometryMode;
  out.coordsAvailable = out.coords.length > 0 && out.atoms.length > 0 && out.coords.length === out.atoms.length;
  out.degenerateCoords = are_atoms_degenerate(out.coords);
  out.usedFallbackCoords = usedFallbackCoords;
  out.ok = out.atoms.length > 0 && out.coordsAvailable;

  if (out.geometryMode === "failed") push_group_warning(out.warningGroups, out.warnings, "geometry", "geometry_failed");

  out.chemistry = assess_smiles_trust({ backend, geometry: out });
  if (backend) backend.chemistry = clone_chemistry_assessment(out.chemistry);

  dlog("SMILES geometry result", {
    ok: out.ok,
    geometryMode: out.geometryMode,
    coordsAvailable: out.coordsAvailable,
    degenerateCoords: out.degenerateCoords,
    usedFallbackCoords: out.usedFallbackCoords,
    warnings: out.warnings
  });

  return out;
}

function finalize_smiles_geometry(input) {
  return build_smiles_geometry_result(input);
}

function assess_smiles_trust(input) {
  let backend = resolve_backend_from_input_for_trust(input);
  let geometry = resolve_geometry_from_input_for_trust(input);
  let ownsBackend = false;

  if (typeof input === "string") {
    backend = build_smiles_backend_result(input);
    ownsBackend = true;
  }

  if (!backend && input && typeof input === "object" && input.mol && input.identity && input.stages) backend = input;
  if (!backend && input && typeof input === "object" && input.__tem_smiles_backend) backend = input.__tem_smiles_backend;
  if (!geometry && input && typeof input === "object" && input.__tem_smiles_geometry) geometry = input.__tem_smiles_geometry;
  if (!geometry && input && typeof input === "object" && typeof input.geometryMode === "string") geometry = input;
  if (!geometry && backend) geometry = build_smiles_geometry_result(backend);

  const out = make_chemistry_assessment_stub();
  const base = backend ? clone_chemistry_assessment(backend.chemistry) : make_chemistry_assessment_stub();
  out.topologyTrusted = !!base.topologyTrusted;
  out.identityTrusted = !!base.identityTrusted;
  out.hydrogenTrusted = !!base.hydrogenTrusted;

  const geometryMode = geometry && typeof geometry.geometryMode === "string" ? geometry.geometryMode : "failed";
  const geometryWarnings = geometry && Array.isArray(geometry.warnings) ? geometry.warnings : [];
  const authoritativeGeometrySource =
    (geometry && geometry.backend && geometry.backend.authoritativeGeometrySource) ||
    (backend && backend.authoritativeGeometrySource) ||
    "fallback";
  const geometryTrusted = !!(
    geometry &&
    geometry.ok &&
    !geometry.degenerateCoords &&
    !geometry.usedFallbackCoords &&
    (
      geometryMode === "native" ||
      geometryMode === "embed3d" ||
      geometryMode === "replacement_2d" ||
      ((geometryMode === "2d_generated") && authoritativeGeometrySource !== "fallback")
    )
  );

  out.geometryTrusted = geometryTrusted;
  out.renderTrusted = !!(out.topologyTrusted && out.identityTrusted && out.hydrogenTrusted && out.geometryTrusted);
  out.strictOk = out.renderTrusted;

  if (base && Array.isArray(base.rejectReasons)) {
    for (let i = 0; i < base.rejectReasons.length; i++) push_reject_reason(out, base.rejectReasons[i]);
  }

  if (!geometryTrusted) {
    if (geometryMode === "2d_generated") push_reject_reason(out, "geometry_mode_2d_generated");
    else if (geometryMode === "fallback_circle") push_reject_reason(out, "geometry_mode_fallback_circle");
    else if (geometryMode === "replacement_2d") push_reject_reason(out, "geometry_mode_replacement_2d_untrusted");
    else push_reject_reason(out, "geometry_failed");
    if (geometry && geometry.degenerateCoords) push_reject_reason(out, "geometry_degenerate");
    if (has_warning(geometryWarnings, "geometry_2d_generation_failed")) push_reject_reason(out, "geometry_2d_generation_failed");
  }

  if (out.strictOk) out.trustLevel = "strict_ok";
  else if (!out.topologyTrusted) out.trustLevel = "rejected";
  else if (!out.identityTrusted && !out.hydrogenTrusted && !out.geometryTrusted) out.trustLevel = "topology_only";
  else if (!out.identityTrusted) out.trustLevel = "identity_degraded";
  else if (!out.hydrogenTrusted) out.trustLevel = "hydrogen_degraded";
  else if (!out.geometryTrusted) out.trustLevel = "geometry_degraded";
  else out.trustLevel = "rejected";

  dlog("SMILES trust", {
    strictOk: out.strictOk,
    trustLevel: out.trustLevel,
    topologyTrusted: out.topologyTrusted,
    identityTrusted: out.identityTrusted,
    hydrogenTrusted: out.hydrogenTrusted,
    geometryTrusted: out.geometryTrusted,
    renderTrusted: out.renderTrusted,
    rejectReasons: out.rejectReasons
  });

  if (ownsBackend && backend && backend.mol) safe_delete(backend.mol);
  return out;
}

function is_strictly_renderable_smiles(input) {
  return !!assess_smiles_trust(input).strictOk;
}

function get_smiles_identity_info(smiles) {
  const result = build_smiles_backend_result(smiles);
  try {
    if (!result || !result.stages || !result.stages.molCreated) return null;
    return clone_identity_info(result.identity);
  } finally {
    if (result && result.mol) safe_delete(result.mol);
  }
}

export const RDKit = {
  setModule(mod) {
    _RDKit = mod;
  },

  replacement_ready,
  replacement_get_status,
  ensure_smiles_backends_ready,
  normalize_symbol,
  symbol_to_Z_safe,
  get_mol_identity_info,
  get_smiles_identity_info,
  build_smiles_backend_result,
  validate_smiles_backend,
  has_usable_smiles_geometry,
  build_smiles_geometry_result,
  finalize_smiles_geometry,
  assess_smiles_trust,
  is_strictly_renderable_smiles,

  ptable_getAtomicNumber(sym) {
    return symbol_to_Z_safe(sym);
  },

  ptable_getRcovalent(Z) {
    try {
      const mod = req();
      if (typeof mod.getRcovalent === "function") {
        const r = mod.getRcovalent(Z);
        if (Number.isFinite(r)) return r;
      }
    } catch (_) {}
    return 0.77;
  },

  raw_get_mol,
  raw_get_molblock,
  raw_get_json,
  raw_get_num_atoms,
  raw_get_num_bonds,
  raw_get_atom_with_idx,
  raw_addHs,
  raw_embed,
  raw_uff_optimize,
  raw_compute2d,
  raw_generate_aligned_coords,

  get_mol(smiles) {
    const mol = raw_get_mol(smiles);
    if (mol) return mol;
    throw new Error("RDKit: get_mol недоступна у цій збірці");
  },

  get_molblock_text,
  get_json_text,
  get_num_atoms_safe,
  get_num_bonds_safe,
  summarize_molblock,
  add_hs_safe,
  get_atoms_bonds_safe,
  get_coordinates_safe,

  get_coordinates(mol) {
    return get_coordinates_safe(mol);
  },

  get_atoms_bonds(mol) {
    return get_atoms_bonds_safe(mol);
  }
};


function wrap_rdkit_ready_with_replacement(promiseLike) {
  if (!promiseLike || typeof promiseLike.then !== "function") return promiseLike;
  if (promiseLike.__tem_replacementWrapped) return promiseLike;
  const wrapped = Promise.resolve(promiseLike).then(async (mod) => {
    try {
      const maybe = replacement_start_loading();
      if (maybe && typeof maybe.then === "function") await maybe;
    } catch (e) {
      _replacementChem.lastError = e;
      dlog("replacement backend unavailable, keeping RDKit path", e);
    }
    return mod;
  });
  wrapped.__tem_replacementWrapped = true;
  return wrapped;
}

const rdkitReadyPromise = (() => {
  if (window.RDKitReady && typeof window.RDKitReady.then === "function") {
    window.RDKitReady = wrap_rdkit_ready_with_replacement(window.RDKitReady);
    return window.RDKitReady;
  }
  if (typeof window.initRDKitModule === "function") {
    const p = window.initRDKitModule({
      locateFile: () => "./rdkit/RDKit_minimal.wasm"
    });
    window.RDKitReady = wrap_rdkit_ready_with_replacement(p);
    return window.RDKitReady;
  }
  throw new Error("RDKit_minimal.js не підключений або підключений занадто пізно.");
})();

rdkitReadyPromise.then((mod) => {
  _RDKit = mod;

  if (mod.Mol && mod.Mol.prototype) {
    const M = mod.Mol.prototype;

    if (typeof M.addHs === "function") M.__tem_native_addHs = M.addHs;
    if (typeof M.add_hs === "function") M.__tem_native_add_hs = M.add_hs;
    if (typeof M.AddHs === "function") M.__tem_native_AddHs = M.AddHs;
    if (typeof M.embedMolecule === "function") M.__tem_native_embedMolecule = M.embedMolecule;
    if (typeof M.EmbedMolecule === "function") M.__tem_native_EmbedMolecule = M.EmbedMolecule;
    if (typeof M.UFFOptimizeMolecule === "function") M.__tem_native_UFFOptimizeMolecule = M.UFFOptimizeMolecule;
    if (typeof M.uffOptimizeMolecule === "function") M.__tem_native_uffOptimizeMolecule = M.uffOptimizeMolecule;
    if (typeof M.getNumAtoms === "function") M.__tem_native_getNumAtoms = M.getNumAtoms;
    if (typeof M.getNumBonds === "function") M.__tem_native_getNumBonds = M.getNumBonds;
    if (typeof M.generate_alignedCoords === "function") M.__tem_native_generate_alignedCoords = M.generate_alignedCoords;
    if (typeof M.generateAlignedCoords === "function") M.__tem_native_generateAlignedCoords = M.generateAlignedCoords;

    M.add_hs = function () {
      return RDKit.add_hs_safe(this).mol;
    };

    M.embed_molecule = function () {
      return RDKit.raw_embed(this);
    };

    M.uff_optimize_molecule = function () {
      return RDKit.raw_uff_optimize(this);
    };

    M.generate_aligned_coords = function () {
      try {
        const fn = raw_method(this, ["__tem_native_generate_alignedCoords", "__tem_native_generateAlignedCoords", "generate_alignedCoords", "generateAlignedCoords"]);
        if (fn) fn();
      } catch (_) {}
    };

    M.get_num_atoms = function () {
      return RDKit.get_num_atoms_safe(this);
    };

    M.get_num_bonds = function () {
      return RDKit.get_num_bonds_safe(this);
    };

    M.get_atom_atomicnum = function (i) {
      return safe_get_atom_atomicnum(this, i);
    };

    M.get_bond_begin = function (k) {
      const rec = safe_get_bond_record(this, k);
      return rec ? rec[0] : 0;
    };

    M.get_bond_end = function (k) {
      const rec = safe_get_bond_record(this, k);
      return rec ? rec[1] : 0;
    };

    M.get_bond_order = function (k) {
      const rec = safe_get_bond_record(this, k);
      return rec ? rec[2] : 1;
    };
  }

  console.log("✅ RDKit WASM ініціалізовано");
}).catch((e) => {
  console.error("❌ Помилка ініціалізації RDKit:", e);
});
