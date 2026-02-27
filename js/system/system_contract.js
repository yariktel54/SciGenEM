// js/system/system_contract.js
// Minimal, non-invasive contract helpers for "System" build results.
// Goals (A2 wave):
//  - keep current public API intact (builders still return [atoms, bonds, title])
//  - provide optional DEV-only validation to catch silent breakages early
//  - avoid changing semantics (especially bonds null vs [])

export function is_dev_mode() {
  // Opt-in only. Any truthy of these enables validation:
  //  - globalThis.TEM_DEV = true
  //  - URL has ?tem_dev=1
  //  - localStorage TEM_DEV="1"
  try {
    if (typeof globalThis !== 'undefined' && globalThis.TEM_DEV === true) return true;
  } catch (_) {}

  try {
    // location is not available in WebWorkers
    if (typeof location !== 'undefined' && location && location.search) {
      const qs = new URLSearchParams(location.search);
      if (qs.get('tem_dev') === '1') return true;
    }
  } catch (_) {}

  try {
    if (typeof localStorage !== 'undefined') {
      if (localStorage.getItem('TEM_DEV') === '1') return true;
    }
  } catch (_) {}

  return false;
}

export function normalize_build_result(result) {
  // Builders are expected to return: [atoms, bonds, title]
  // We only normalize safe things (types), NOT semantics.
  if (!Array.isArray(result) || result.length < 3) {
    return { atoms: [], bonds: null, title: 'ERROR' };
  }

  const atoms = Array.isArray(result[0]) ? result[0] : [];
  const bonds = result[1]; // keep as-is (null vs [] matters!)
  const title = (result[2] == null) ? '' : String(result[2]);

  return { atoms, bonds, title };
}

export function validate_system(sys, stage = 'build') {
  // DEV-only validation. Throws on hard errors.
  if (!is_dev_mode()) return true;

  const errors = [];
  if (!sys || typeof sys !== 'object') errors.push('system is not an object');

  const atoms = sys?.atoms;
  const bonds = sys?.bonds;

  if (!Array.isArray(atoms)) {
    errors.push('atoms is not an array');
  } else {
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      if (!a || typeof a !== 'object') {
        errors.push(`atom[${i}] is not an object`);
        break;
      }
      const x = a.x;
      const y = a.y;
      const z = ('z' in a) ? a.z : 0;
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
        errors.push(`atom[${i}] has non-finite coords (x,y,z)`);
        break;
      }
      if (!('Z' in a) || !Number.isFinite(a.Z)) {
        // Not a hard error for some pipelines, but very useful to know.
        errors.push(`atom[${i}] has missing/invalid Z`);
        break;
      }
    }
  }

  if (!(bonds === null || Array.isArray(bonds))) {
    errors.push('bonds must be null or an array');
  }

  if (Array.isArray(bonds) && Array.isArray(atoms)) {
    const n = atoms.length;
    for (let k = 0; k < bonds.length; k++) {
      const b = bonds[k];
      if (!Array.isArray(b) || b.length < 2) {
        errors.push(`bond[${k}] is not an [i,j,...] array`);
        break;
      }
      const i = b[0] | 0;
      const j = b[1] | 0;
      if (i < 0 || j < 0 || i >= n || j >= n) {
        errors.push(`bond[${k}] indices out of range (n=${n})`);
        break;
      }
      if (i === j) {
        errors.push(`bond[${k}] has i===j`);
        break;
      }
    }
  }

  if (errors.length) {
    throw new Error(`[TEM][${stage}] System validation failed: ${errors[0]}`);
  }

  return true;
}
