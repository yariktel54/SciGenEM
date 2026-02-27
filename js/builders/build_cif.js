// js/builders/build_cif.js
// CIF system builder (tiling + bonds semantics) used by phases.js.
// Safe refactor: logic matches the previous implementation in phases.js.

import { load_cif } from '../io/cif_io.js';
import { guess_bonds_by_distance } from '../bonds/bonds.js';
import { cell_vectors, frac_to_cart, tile_shift } from '../system/lattice.js';

function normalize_size(size) {
  if (!size) return [1, 1, 1];
  if (size === 'auto') return [1, 1, 1];
  if (!Array.isArray(size) || size.length !== 3) return [1, 1, 1];
  const nx = (size[0] | 0) || 1;
  const ny = (size[1] | 0) || 1;
  const nz = (size[2] | 0) || 1;
  return [nx, ny, nz];
}

export async function build_from_cif(spec) {
  const data = await load_cif(spec.path);

  // Cell params are required for fractional conversion and for tiling shifts.
  // data.cell is expected to be [a, b, c, alpha, beta, gamma]
  const cell = data?.cell;
  if (!cell || cell.length < 6) {
    // If CIF has no cell, we cannot do fractional->cart or proper tiling.
    // Still try to return any available cart coords as-is.
    const atoms_cart = Array.isArray(data?.atoms_cart) ? data.atoms_cart : [];
    const atoms_fallback = atoms_cart.map(at => ({ Z: at.Z, x: at.x, y: at.y, z: at.z }));
    const bonds_fallback = (data && Array.isArray(data.bonds)) ? data.bonds : null;
    return [atoms_fallback, bonds_fallback, 'CIF'];
  }

  const a = cell[0], b = cell[1], c = cell[2], alpha = cell[3], beta = cell[4], gamma = cell[5];
  const basis = cell_vectors(a, b, c, alpha, beta, gamma);
  const vx = basis[0], vy = basis[1], vz = basis[2];

  const size = normalize_size(spec.size);
  const nx = size[0], ny = size[1], nz = size[2];

  // Build one-cell cart coords.
  const unit = [];

  const frac = Array.isArray(data?.atoms_frac) ? data.atoms_frac : [];
  const cart = Array.isArray(data?.atoms_cart) ? data.atoms_cart : [];

  if (frac.length > 0) {
    for (const at of frac) {
      const xyz = frac_to_cart([at.fx, at.fy, at.fz], vx, vy, vz);
      unit.push({ Z: at.Z, x: xyz[0], y: xyz[1], z: xyz[2] });
    }
  } else if (cart.length > 0) {
    for (const at of cart) {
      unit.push({ Z: at.Z, x: at.x, y: at.y, z: at.z });
    }
  } else {
    // Nothing parsed
    const bonds_empty = (data && Array.isArray(data.bonds)) ? data.bonds : null;
    return [[], bonds_empty, 'CIF (no atoms)'];
  }

  // Tile if needed
  const atoms = [];
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      for (let iz = 0; iz < nz; iz++) {
        const d = tile_shift(ix, iy, iz, vx, vy, vz);
        const dx = d[0], dy = d[1], dz = d[2];
        for (const at of unit) {
          atoms.push({ Z: at.Z, x: at.x + dx, y: at.y + dy, z: at.z + dz });
        }
      }
    }
  }

  // --- CIF metadata (cell is useful for debugging/scale bars etc) ---
  atoms._cell = [a, b, c, alpha, beta, gamma];

  // IMPORTANT:
  //  - we already explicitly tile atoms (nx,ny,nz), so we MUST NOT apply PBC minimum-image on the tiled supercell
  //    (that creates wrap-around "ghost" short distances and extra bonds).
  //  - for nx=ny=nz=1 we also keep PBC OFF to avoid "long bonds" drawn across the cell boundary;
  //    if periodic connectivity is needed, tile to >1 and we will catch neighbours by plain distance.
  atoms._pbc = false;
  delete atoms._cell_vectors;

  const unitN = unit.length;

  // 1) bonds from CIF tables (if present) are explicit → we must not guess.
  //    If we tiled the system, replicate intra-cell bonds to each tile.
  let bonds = (data && Array.isArray(data.bonds)) ? data.bonds : null;
  if (Array.isArray(bonds)) {
    if (nx * ny * nz > 1 && unitN > 0 && bonds.length > 0) {
      const tiled = [];
      let tileIndex = 0;
      for (let ix = 0; ix < nx; ix++) {
        for (let iy = 0; iy < ny; iy++) {
          for (let iz = 0; iz < nz; iz++) {
            const off = tileIndex * unitN;
            for (const bnd of bonds) {
              tiled.push([(bnd[0] | 0) + off, (bnd[1] | 0) + off, bnd[2] ?? 1]);
            }
            tileIndex++;
          }
        }
      }
      bonds = tiled;
    }
  }

  // 2) If CIF has no bonds table → guess once here (bonds.js is the single source of truth).
  if (bonds === null) {
    bonds = guess_bonds_by_distance(atoms, 1.22, 2.1);
    if (!Array.isArray(bonds)) bonds = [];
  }

  const title = (nx === 1 && ny === 1 && nz === 1) ? 'CIF' : `CIF size=${nx}x${ny}x${nz}`;
  return [atoms, bonds, title];
}
