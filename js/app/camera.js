// js/app/camera.js
// Minimal camera / view math for TEM simulator.

export function make_view_state(opts) {
  opts = opts || {};
  const img_size = Array.isArray(opts.img_size) ? opts.img_size : [400, 400];
  const H = img_size[0] | 0;
  const W = img_size[1] | 0;

  let ang = Number(opts.angstroms_per_pixel);
  if (!Number.isFinite(ang) || ang <= 0) ang = 0.1;

  const pan_px = Array.isArray(opts.pan_px) ? opts.pan_px : [0, 0];
  let panX = Number(pan_px[0]);
  let panY = Number(pan_px[1]);
  if (!Number.isFinite(panX)) panX = 0;
  if (!Number.isFinite(panY)) panY = 0;

  let rotZ = Number(opts.rotZ_rad);
  if (!Number.isFinite(rotZ)) rotZ = 0;

  const center_A = Array.isArray(opts.center_A) ? opts.center_A : [0, 0, 0];
  let cxA = Number(center_A[0]);
  let cyA = Number(center_A[1]);
  let czA = Number(center_A[2]);
  if (!Number.isFinite(cxA)) cxA = 0;
  if (!Number.isFinite(cyA)) cyA = 0;
  if (!Number.isFinite(czA)) czA = 0;

  return {
    img_size: [H, W],
    angstroms_per_pixel: ang,
    pan_px: [panX, panY],
    rotZ_rad: rotZ,
    center_A: [cxA, cyA, czA],
    center_mode: opts.center_mode || 'bbox'
  };
}

function _screen_center_px(viewState) {
  const H = viewState.img_size[0];
  const W = viewState.img_size[1];
  return [(W - 1) * 0.5, (H - 1) * 0.5];
}

export function unproject_xy(viewState, sx, sy) {
  const sc = _screen_center_px(viewState);
  const cx = sc[0] + viewState.pan_px[0];
  const cy = sc[1] + viewState.pan_px[1];
  const dx = (sx - cx) * viewState.angstroms_per_pixel;
  const dy = (sy - cy) * viewState.angstroms_per_pixel;

  const a = -(viewState.rotZ_rad || 0);
  const ca = Math.cos(a);
  const sa = Math.sin(a);
  const xr = dx * ca - dy * sa;
  const yr = dx * sa + dy * ca;

  return [xr + viewState.center_A[0], yr + viewState.center_A[1]];
}

// World-space AABB of screen rectangle in Å (XY only) with optional padding in px.
export function world_aabb_xy(viewState, pad_px) {
  pad_px = Number(pad_px);
  if (!Number.isFinite(pad_px)) pad_px = 0;

  const H = viewState.img_size[0];
  const W = viewState.img_size[1];

  const c0 = unproject_xy(viewState, -pad_px, -pad_px);
  const c1 = unproject_xy(viewState, W - 1 + pad_px, -pad_px);
  const c2 = unproject_xy(viewState, -pad_px, H - 1 + pad_px);
  const c3 = unproject_xy(viewState, W - 1 + pad_px, H - 1 + pad_px);

  const xs = [c0[0], c1[0], c2[0], c3[0]];
  const ys = [c0[1], c1[1], c2[1], c3[1]];

  return {
    minx: Math.min.apply(null, xs),
    maxx: Math.max.apply(null, xs),
    miny: Math.min.apply(null, ys),
    maxy: Math.max.apply(null, ys)
  };
}
