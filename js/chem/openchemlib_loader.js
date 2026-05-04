// js/chem/openchemlib_loader.js
// Loader for ESM build of OpenChemLib.

const OPEN_CHEM_LIB_VENDOR_RELATIVE_PATH = "./vendor/openchemlib.js";
const OPEN_CHEM_LIB_VENDOR_URL = new URL(
  OPEN_CHEM_LIB_VENDOR_RELATIVE_PATH,
  import.meta.url,
).href;

let _module = null;
let _ready = null;
let _lastError = null;
let _source = "none";

function normalizeOpenChemLibModule(mod) {
  const candidates = [mod, mod && mod.default, mod && mod.OCL];

  for (let i = 0; i < candidates.length; i++) {
    const c = candidates[i];
    if (c && c.Molecule && typeof c.Molecule.fromSmiles === "function") {
      return c;
    }
  }

  try {
    if (typeof window !== "undefined" && window && window.OCL) {
      const c = window.OCL;
      if (c && c.Molecule && typeof c.Molecule.fromSmiles === "function") {
        return c;
      }
    }
  } catch (_) {}

  return null;
}

export async function ensureOpenChemLibLoaded() {
  if (_module) return _module;
  if (_ready) return _ready;

  _ready = import(/* @vite-ignore */ OPEN_CHEM_LIB_VENDOR_URL)
    .then((mod) => {
      const normalized = normalizeOpenChemLibModule(mod);
      if (!normalized) {
        throw new Error("openchemlib_vendor_invalid_module");
      }
      _module = normalized;
      _source = OPEN_CHEM_LIB_VENDOR_URL;
      _lastError = null;
      return normalized;
    })
    .catch((e) => {
      _lastError = e;
      _ready = null;
      throw e;
    });

  return _ready;
}

export function getOpenChemLibModule() {
  if (_module) return _module;
  return null;
}

export function getOpenChemLibLoaderStatus() {
  return {
    ready: !!_module,
    source: _source || OPEN_CHEM_LIB_VENDOR_URL,
    error: _lastError ? String(_lastError.message || _lastError) : "",
  };
}
