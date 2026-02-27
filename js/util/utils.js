// js/util/utils.js
export function normalize_and_scale(coords, size) {
    // (не використовується напряму в цій JS версії, але збережено для сумісності)
    const n = coords.length;
    const mean = coords.reduce((acc, [x, y]) => [acc[0] + x, acc[1] + y], [0, 0]).map(s => s / n);
    const centered = coords.map(([x, y]) => [x - mean[0], y - mean[1]]);
    let maxRange = 0;
    for (const [x, y] of centered) maxRange = Math.max(maxRange, Math.abs(x), Math.abs(y));
    const scale = (size * 0.4) / Math.max(maxRange, 1e-9);
    return centered.map(([x, y]) => [x * scale + size / 2, y * scale + size / 2]);
}
