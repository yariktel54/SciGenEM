// js/system/lattice.js
// Crystallographic helpers (cell vectors, fractional/cartesian conversion, tiling shifts).
// Pure math; no IO and no domain logic.

export function cell_vectors(a, b, c, alpha = 90, beta = 90, gamma = 90) {
  const ar = alpha * Math.PI / 180;
  const br = beta * Math.PI / 180;
  const gr = gamma * Math.PI / 180;

  // Standard crystallographic basis vectors
  const vx = [a, 0, 0];
  const vy = [b * Math.cos(gr), b * Math.sin(gr), 0];

  const cx = c * Math.cos(br);
  const cy = c * (Math.cos(ar) - Math.cos(br) * Math.cos(gr)) / Math.max(1e-12, Math.sin(gr));
  const cz2 = c * c - cx * cx - cy * cy;
  const cz = Math.sqrt(Math.max(0, cz2));
  const vz = [cx, cy, cz];

  return [vx, vy, vz];
}

export function frac_to_cart(f, vx, vy, vz) {
  const fx = f[0], fy = f[1], fz = f[2];
  return [
    fx * vx[0] + fy * vy[0] + fz * vz[0],
    fx * vx[1] + fy * vy[1] + fz * vz[1],
    fx * vx[2] + fy * vy[2] + fz * vz[2]
  ];
}

export function tile_shift(ix, iy, iz, vx, vy, vz) {
  return [
    ix * vx[0] + iy * vy[0] + iz * vz[0],
    ix * vx[1] + iy * vy[1] + iz * vz[1],
    ix * vx[2] + iy * vy[2] + iz * vz[2]
  ];
}
