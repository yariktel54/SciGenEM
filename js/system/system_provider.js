// js/system/system_provider.js
// Provider layer for "true lazy" view-dependent subsets.
//
// V1 scope:
//  - Stop CIF eager tiling by using PeriodicProvider.
//  - Keep renderer compatible: getView() returns atomsView as SMALL array of {Z,x,y,z} objects.

import { world_aabb_xy } from '../app/camera.js';
import { guess_bonds_by_distance } from '../bonds/bonds_guess.js';
import { cell_vectors, frac_to_cart, tile_shift } from '../system/lattice.js';

function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

function build_cart_to_frac(vx, vy, vz) {
  const bxc = cross(vy, vz);
  const cxa = cross(vz, vx);
  const axb = cross(vx, vy);
  const det = dot(vx, bxc);
  if (!Number.isFinite(det) || Math.abs(det) < 1e-10) return null;

  const inv = [
    [bxc[0] / det, bxc[1] / det, bxc[2] / det],
    [cxa[0] / det, cxa[1] / det, cxa[2] / det],
    [axb[0] / det, axb[1] / det, axb[2] / det]
  ];
  return function cart_to_frac(v) {
    return [
      inv[0][0] * v[0] + inv[0][1] * v[1] + inv[0][2] * v[2],
      inv[1][0] * v[0] + inv[1][1] * v[1] + inv[1][2] * v[2],
      inv[2][0] * v[0] + inv[2][1] * v[1] + inv[2][2] * v[2]
    ];
  };
}

function bbox_from_atoms(atoms) {
  let xmin = Infinity, ymin = Infinity, zmin = Infinity;
  let xmax = -Infinity, ymax = -Infinity, zmax = -Infinity;
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    const x = a.x, y = a.y, z = a.z || 0;
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) continue;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }
  if (!Number.isFinite(xmin)) xmin = ymin = zmin = xmax = ymax = zmax = 0;
  return { xmin, xmax, ymin, ymax, zmin, zmax };
}

function center_from_bbox(bb) {
  return [0.5 * (bb.xmin + bb.xmax), 0.5 * (bb.ymin + bb.ymax), 0.5 * (bb.zmin + bb.zmax)];
}

function clamp_int(v, lo, hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

function tiles_allowed_range(size) {
  if (size === 'auto') {
    return { finite: false, ix0: -Infinity, ix1: Infinity, iy0: -Infinity, iy1: Infinity, iz0: -Infinity, iz1: Infinity, nx: null, ny: null, nz: null };
  }
  const nx = Math.max(1, size[0] | 0);
  const ny = Math.max(1, size[1] | 0);
  const nz = Math.max(1, size[2] | 0);
  return { finite: true, ix0: 0, ix1: nx - 1, iy0: 0, iy1: ny - 1, iz0: 0, iz1: nz - 1, nx, ny, nz };
}

// -------- Bonds fallback (single implementation) --------
// Semantics (strict):
//  - bonds === null  -> unknown -> guess allowed (as fallback) BUT only on the *view subset*.
//  - bonds is array  -> explicit topology (even []) -> guess forbidden.
//
// We cache guessed bonds per provider+view to avoid recomputation on pan/zoom/rot.
function _q(v, s) { return Math.round((Number(v) || 0) * s); }
function _qInt(v) { return (Number.isFinite(v) ? (v | 0) : 0); }

function _make_view_key(sysHash, cam, extra) {
  cam = cam || {};
  const img = Array.isArray(cam.img_size) ? cam.img_size : [0, 0]; // [H,W]
  const H = _qInt(img[0]);
  const W = _qInt(img[1]);

  const ap = _q(cam.angstroms_per_pixel, 10000);
  const pan = Array.isArray(cam.pan_px) ? cam.pan_px : [0, 0];
  const panX = _qInt(pan[0]);
  const panY = _qInt(pan[1]);

  const rz = _q(cam.rotZ_rad, 10000);
  const rx = _q(cam.rotX_rad, 10000);
  const ry = _q(cam.rotY_rad, 10000);

  const c = Array.isArray(cam.center_A) ? cam.center_A : [0, 0, 0];
  const cx = _q(c[0], 1000);
  const cy = _q(c[1], 1000);
  const cz = _q(c[2], 1000);

  let k = `${sysHash}|${W}x${H}|ap=${ap}|pan=${panX},${panY}|rz=${rz}|rx=${rx}|ry=${ry}|c=${cx},${cy},${cz}`;
  if (extra) k += '|' + extra;
  return k;
}

function _cache_put_lru(map, key, val, maxSize) {
  map.set(key, val);
  if (map.size <= maxSize) return;
  const firstKey = map.keys().next().value;
  if (firstKey !== undefined) map.delete(firstKey);
}

function _mark_guessed(arr) {
  try {
    if (arr && typeof arr === 'object') arr.__tem_guessed = true;
  } catch (_) {}
  return arr;
}



// -------- StaticProvider --------
export class StaticProvider {
  constructor(params) {
    params = params || {};
    this._atoms = Array.isArray(params.atoms) ? params.atoms : [];
    this._bonds = (params.bonds === null) ? null : (Array.isArray(params.bonds) ? params.bonds : []);
    this._title = params.title || '';
    this._meta = params.meta || {};
    this._bbox = bbox_from_atoms(this._atoms);
    this._center_A = center_from_bbox(this._bbox);

    // guess cache (cleared automatically on rebuild because provider instance is replaced)
    this._guessCache = new Map();
    this._sysHash = 'S|' + this._atoms.length + '|' + _q(this._bbox.xmin, 100) + ',' + _q(this._bbox.ymin, 100) + ',' + _q(this._bbox.zmin, 100) + '|' + _q(this._bbox.xmax, 100) + ',' + _q(this._bbox.ymax, 100) + ',' + _q(this._bbox.zmax, 100);
  }

  getMeta() {
    return {
      title: this._title,
      kind: 'molecule',
      hasBonds: this._bonds !== null,
      bondsKnown: this._bonds !== null,
      atomCount: this._atoms.length,
      center_A: this._center_A,
      bbox: this._bbox,
      meta: this._meta
    };
  }

  getView(viewState, opts) {
    opts = opts || {};
    const usedCamera = Object.assign({}, viewState || {});
    usedCamera.center_A = this._center_A;

    const needBonds = (opts.needBonds !== false);
    let bondsView = this._bonds;

    // If bonds are unknown (null) -> guess allowed as a late fallback (view-dependent).
    if (this._bonds === null) {
      if (!needBonds) {
        bondsView = null;
      } else {
        const key = _make_view_key(this._sysHash, usedCamera, 'static');
        const cached = this._guessCache.get(key);
        if (cached !== undefined) {
          // refresh LRU
          this._guessCache.delete(key);
          this._guessCache.set(key, cached);
          bondsView = cached;
        } else {
          // Guess only for atoms likely in-frame to avoid O(N^2) on huge systems.
          const aabb = world_aabb_xy(usedCamera, 96);
          const atoms = this._atoms;

          const idxs = [];
          for (let i = 0; i < atoms.length; i++) {
            const a = atoms[i];
            const x = a.x, y = a.y;
            if (x < aabb.minx || x > aabb.maxx || y < aabb.miny || y > aabb.maxy) continue;
            idxs.push(i);
          }

          const MAX_GUESS_ATOMS = 8000;
          if (idxs.length >= 2) {
            // If there are too many atoms in-frame, downsample deterministically.
            let pick = idxs;
            if (pick.length > MAX_GUESS_ATOMS) {
              const step = Math.ceil(pick.length / MAX_GUESS_ATOMS);
              const sampled = [];
              for (let t = 0; t < pick.length; t += step) sampled.push(pick[t]);
              pick = sampled;
            }

            const sub = new Array(pick.length);
            for (let k = 0; k < pick.length; k++) sub[k] = atoms[pick[k]];

            const guessedLocal = guess_bonds_by_distance(sub, { periodic: false, tolA: 0.45, maxNeighborsPerAtom: 12 });
            const guessed = new Array(guessedLocal.length);
            for (let k = 0; k < guessedLocal.length; k++) {
              const b = guessedLocal[k];
              guessed[k] = [pick[b[0]], pick[b[1]], b[2] || 1];
            }
            bondsView = _mark_guessed(guessed);
          } else {
            // Skip guessing for this view (empty).
            bondsView = [];
          }

          _cache_put_lru(this._guessCache, key, bondsView, 64);
        }
      }
    }

    return { atomsView: this._atoms, bondsView, title: this._title, usedCamera };
  }
}

// -------- PeriodicProvider (CIF) --------
export class PeriodicProvider {
  constructor(params) {
    params = params || {};
    this._title = params.title || '';
    this._cell = params.cell || null; // [a,b,c,alpha,beta,gamma]
    this._size = (params.size === 'auto') ? 'auto' : (Array.isArray(params.size) ? params.size : [1, 1, 1]);
    this._baseBonds = (params.bonds === null) ? null : (Array.isArray(params.bonds) ? params.bonds : []);

    const a = this._cell ? this._cell[0] : 1;
    const b = this._cell ? this._cell[1] : 1;
    const c = this._cell ? this._cell[2] : 1;
    const alpha = this._cell ? this._cell[3] : 90;
    const beta = this._cell ? this._cell[4] : 90;
    const gamma = this._cell ? this._cell[5] : 90;

    const vv = cell_vectors(a, b, c, alpha, beta, gamma);
    this._vx = vv[0]; this._vy = vv[1]; this._vz = vv[2];
    this._cart_to_frac = build_cart_to_frac(this._vx, this._vy, this._vz);

    const atoms_cart_in = Array.isArray(params.atoms_cart) ? params.atoms_cart : [];
    const atoms_frac_in = Array.isArray(params.atoms_frac) ? params.atoms_frac : [];

    // Canonical base atoms for periodic systems:
    //  - prefer fractional coordinates when available (stable + matches CIF)
    //  - NEVER merge/concatenate cart+frac (avoids duplicates / ghost atoms)
    const useFrac = (atoms_frac_in.length > 0);
    const n = useFrac ? atoms_frac_in.length : atoms_cart_in.length;
    this._nBase = n;
    this._atomMode = useFrac ? 'frac' : 'cart';

    // DEV: catch accidental dual arrays with mismatched lengths early
    if (globalThis && globalThis.__TEM_DEV_MODE) {
      try {
        if (atoms_cart_in.length && atoms_frac_in.length && atoms_cart_in.length !== atoms_frac_in.length) {
          console.debug(`[TEM] WARN periodic atoms mismatch: atoms_cart=${atoms_cart_in.length}, atoms_frac=${atoms_frac_in.length} (using ${this._atomMode})`);
        }
      } catch (_) {}
    }

    // base cell in SoA (typed arrays) — cheap to keep, fast to tile
    this._Z = new Uint16Array(n);
    this._x = new Float32Array(n);
    this._y = new Float32Array(n);
    this._z = new Float32Array(n);

    if (useFrac) {
      for (let i = 0; i < n; i++) {
        const af = atoms_frac_in[i];
        this._Z[i] = (af && af.Z) ? (af.Z | 0) : 6;
        const p = frac_to_cart([af.fx || 0, af.fy || 0, af.fz || 0], this._vx, this._vy, this._vz);
        this._x[i] = p[0]; this._y[i] = p[1]; this._z[i] = p[2];
      }
    } else {
      for (let i = 0; i < n; i++) {
        const a0 = atoms_cart_in[i];
        this._Z[i] = (a0 && a0.Z) ? (a0.Z | 0) : 6;
        this._x[i] = a0 ? (a0.x || 0) : 0;
        this._y[i] = a0 ? (a0.y || 0) : 0;
        this._z[i] = a0 ? (a0.z || 0) : 0;
      }
    }

    // bbox/center of base unit cell
    const tmp = new Array(n);
    for (let i = 0; i < n; i++) tmp[i] = { x: this._x[i], y: this._y[i], z: this._z[i] };
    this._baseBBox = bbox_from_atoms(tmp);
    this._baseCenter = center_from_bbox(this._baseBBox);

    // stable camera anchor: center of tiling domain (finite) or base-center (auto)
    const lim = tiles_allowed_range(this._size);
    if (lim.finite) {
      const sx = (lim.nx - 1) * 0.5;
      const sy = (lim.ny - 1) * 0.5;
      const sz = (lim.nz - 1) * 0.5;
      const sh = tile_shift(sx, sy, sz, this._vx, this._vy, this._vz);
      this._center_A = [this._baseCenter[0] + sh[0], this._baseCenter[1] + sh[1], this._baseCenter[2] + sh[2]];
    } else {
      this._center_A = [this._baseCenter[0], this._baseCenter[1], this._baseCenter[2]];
    }

    this._didLogView = false;

    // guess cache for view-dependent bond perception (only when base bonds are unknown)
    this._guessCache = new Map();
    this._sysHash = 'P|' + this._nBase + '|' + (this._cell ? (this._cell.map(v => _q(v, 100)).join(',')) : 'nocell') + '|' + _q(this._baseCenter[0], 100) + ',' + _q(this._baseCenter[1], 100) + ',' + _q(this._baseCenter[2], 100);
  }

  getMeta() {
    return {
      title: this._title,
      kind: 'periodic',
      hasBonds: this._baseBonds !== null,
      bondsKnown: this._baseBonds !== null,
      atomCountBase: this._nBase,
      cell: this._cell,
      size: this._size,
      center_A: this._center_A,
      baseBBox: this._baseBBox
    };
  }

  getView(viewState, opts) {
    opts = opts || {};
    const usedCamera = Object.assign({}, viewState || {});
    usedCamera.center_A = this._center_A;

    // screen→world bbox → approximate tile range in fractional coords
    const pad_px = 96;
    const aabb = world_aabb_xy(usedCamera, pad_px);
    const cart_to_frac = this._cart_to_frac;

    let ixMin = 0, ixMax = 0, iyMin = 0, iyMax = 0, izMin = 0, izMax = 0;

    if (typeof cart_to_frac === 'function') {
      const corners = [
        [aabb.minx, aabb.miny, this._center_A[2]],
        [aabb.maxx, aabb.miny, this._center_A[2]],
        [aabb.minx, aabb.maxy, this._center_A[2]],
        [aabb.maxx, aabb.maxy, this._center_A[2]]
      ];
      let fxMin = Infinity, fxMax = -Infinity, fyMin = Infinity, fyMax = -Infinity, fzMin = Infinity, fzMax = -Infinity;
      for (let i = 0; i < corners.length; i++) {
        const f = cart_to_frac(corners[i]);
        if (!f || f.length !== 3) continue;
        if (f[0] < fxMin) fxMin = f[0];
        if (f[0] > fxMax) fxMax = f[0];
        if (f[1] < fyMin) fyMin = f[1];
        if (f[1] > fyMax) fyMax = f[1];
        if (f[2] < fzMin) fzMin = f[2];
        if (f[2] > fzMax) fzMax = f[2];
      }
      if (Number.isFinite(fxMin)) {
        ixMin = Math.floor(fxMin) - 1; ixMax = Math.ceil(fxMax) + 1;
        iyMin = Math.floor(fyMin) - 1; iyMax = Math.ceil(fyMax) + 1;
        izMin = Math.floor(fzMin) - 1; izMax = Math.ceil(fzMax) + 1;
      }
    }

    // clamp to user-specified finite size
    const lim = tiles_allowed_range(this._size);
    if (lim.finite) {
      ixMin = clamp_int(ixMin, lim.ix0, lim.ix1); ixMax = clamp_int(ixMax, lim.ix0, lim.ix1);
      iyMin = clamp_int(iyMin, lim.iy0, lim.iy1); iyMax = clamp_int(iyMax, lim.iy0, lim.iy1);
      izMin = clamp_int(izMin, lim.iz0, lim.iz1); izMax = clamp_int(izMax, lim.iz0, lim.iz1);
    }

    // hard cap: prevent "zoomed out" exploding tiles
    const MAX_TILES_TOTAL = 4096;
    let tilesTotal = (ixMax - ixMin + 1) * (iyMax - iyMin + 1) * (izMax - izMin + 1);
    if (tilesTotal > MAX_TILES_TOTAL) {
      let cx = 0, cy = 0, cz = 0;
      if (typeof cart_to_frac === 'function') {
        const fc = cart_to_frac(this._center_A);
        cx = Math.floor(fc[0]); cy = Math.floor(fc[1]); cz = Math.floor(fc[2]);
      }
      const span = Math.floor(Math.cbrt(MAX_TILES_TOTAL));
      const h = Math.max(1, Math.floor(span / 2));
      ixMin = cx - h; ixMax = cx + h;
      iyMin = cy - h; iyMax = cy + h;
      izMin = cz - h; izMax = cz + h;
      if (lim.finite) {
        ixMin = clamp_int(ixMin, lim.ix0, lim.ix1); ixMax = clamp_int(ixMax, lim.ix0, lim.ix1);
        iyMin = clamp_int(iyMin, lim.iy0, lim.iy1); iyMax = clamp_int(iyMax, lim.iy0, lim.iy1);
        izMin = clamp_int(izMin, lim.iz0, lim.iz1); izMax = clamp_int(izMax, lim.iz0, lim.iz1);
      }
    }

    // coarse tile bbox filter in world XY
    const tiles = [];
    const bb = this._baseBBox;
    for (let iz = izMin; iz <= izMax; iz++) {
      for (let iy = iyMin; iy <= iyMax; iy++) {
        for (let ix = ixMin; ix <= ixMax; ix++) {
          const sh = tile_shift(ix, iy, iz, this._vx, this._vy, this._vz);
          const tminx = bb.xmin + sh[0], tmaxx = bb.xmax + sh[0];
          const tminy = bb.ymin + sh[1], tmaxy = bb.ymax + sh[1];
          if (tmaxx < aabb.minx || tminx > aabb.maxx || tmaxy < aabb.miny || tminy > aabb.maxy) continue;
          tiles.push({ ix, iy, iz, sh });
        }
      }
    }

    // materialize atoms for visible tiles only (renderer-compatible objects)
    const atomsView = [];
    const nBase = this._nBase;
    for (let t = 0; t < tiles.length; t++) {
      const sh = tiles[t].sh;
      for (let i = 0; i < nBase; i++) {
        atomsView.push({ Z: this._Z[i], x: this._x[i] + sh[0], y: this._y[i] + sh[1], z: this._z[i] + sh[2] });
      }
    }

    // bonds policy (strict semantics):
    //  - baseBonds === null -> unknown -> guess allowed as view-dependent fallback (ONLY if needBonds)
    //  - baseBonds is array (even []) -> explicit topology -> guess forbidden
    let bondsView = this._baseBonds;
    const needBonds = (opts.needBonds !== false);

    if (this._baseBonds === null) {
      if (!needBonds) {
        bondsView = null;
      } else {
        const extraKey = 'tiles=' + ixMin + ',' + ixMax + ',' + iyMin + ',' + iyMax + ',' + izMin + ',' + izMax;
        const key = _make_view_key(this._sysHash, usedCamera, extraKey);
        const cached = this._guessCache.get(key);
        if (cached !== undefined) {
          // refresh LRU
          this._guessCache.delete(key);
          this._guessCache.set(key, cached);
          bondsView = cached;
        } else {
          // Guess bonds only on atoms likely in-frame (avoids O(N^2) on large tiled views).
          const aabb2 = world_aabb_xy(usedCamera, 64);
          const idxs = [];
          for (let i = 0; i < atomsView.length; i++) {
            const a = atomsView[i];
            const x = a.x, y = a.y;
            if (x < aabb2.minx || x > aabb2.maxx || y < aabb2.miny || y > aabb2.maxy) continue;
            idxs.push(i);
          }

          const MAX_GUESS_ATOMS = 8000;
          if (idxs.length >= 2) {
            let pick = idxs;
            if (pick.length > MAX_GUESS_ATOMS) {
              const step = Math.ceil(pick.length / MAX_GUESS_ATOMS);
              const sampled = [];
              for (let t = 0; t < pick.length; t += step) sampled.push(pick[t]);
              pick = sampled;
            }

            const sub = new Array(pick.length);
            for (let k = 0; k < pick.length; k++) sub[k] = atomsView[pick[k]];

            const guessedLocal = guess_bonds_by_distance(sub, { periodic: false, tolA: 0.45, maxNeighborsPerAtom: 16 });
            const guessed = new Array(guessedLocal.length);
            for (let k = 0; k < guessedLocal.length; k++) {
              const b = guessedLocal[k];
              guessed[k] = [pick[b[0]], pick[b[1]], b[2] || 1];
            }
            bondsView = _mark_guessed(guessed);
          } else {
            bondsView = [];
          }
          _cache_put_lru(this._guessCache, key, bondsView, 64);
        }
      }
    } else if (Array.isArray(this._baseBonds) && this._baseBonds.length === 0) {
      bondsView = [];
    } else {
      // Explicit bonds: replicate per tile (can be skipped for huge views).
      const MAX_BOND_ATOMS = 60000;
      if (!needBonds || atomsView.length > MAX_BOND_ATOMS) {
        bondsView = []; // explicit-but-skipped => forbids guessing
      } else {
        const out = [];
        for (let t = 0; t < tiles.length; t++) {
          const off = t * nBase;
          for (let k = 0; k < this._baseBonds.length; k++) {
            const b0 = this._baseBonds[k];
            out.push([b0[0] + off, b0[1] + off, b0[2] || 1]);
          }
        }
        bondsView = out;
      }
    }

    usedCamera._tiles_in_view = tiles.length;
    usedCamera._atoms_in_view = atomsView.length;
    // DEV: single-line view summary (no spam)
    if (!this._didLogView && globalThis && globalThis.__TEM_DEV_MODE) {
      try {
        this._didLogView = true;
        const nb = (bondsView === null) ? 'null' : String(Array.isArray(bondsView) ? bondsView.length : 0);
        const c = this._cell;
        const cellStr = Array.isArray(c) ? c.slice(0,6).map(v => (Number.isFinite(+v) ? (+v).toFixed(6) : String(v))).join(',') : 'null';
        const ca = usedCamera && Array.isArray(usedCamera.center_A) ? usedCamera.center_A : this._center_A;
        const centerStr = Array.isArray(ca) ? ca.map(v => (Number.isFinite(+v) ? (+v).toFixed(6) : String(v))).join(',') : 'null';
        const tilesN = usedCamera && Number.isFinite(usedCamera._tiles_in_view) ? usedCamera._tiles_in_view : null;
        console.debug(`[TEM] view periodic(${this._title || 'periodic'}): atomsView=${atomsView.length}; bondsView=${nb}; cell=${cellStr}; center_A=${centerStr}` + (tilesN != null ? `; tiles=${tilesN}` : '') + `; atomMode=${this._atomMode}`);
      } catch (_) {}
    }

    return { atomsView, bondsView, title: this._title, usedCamera };
  }
}
