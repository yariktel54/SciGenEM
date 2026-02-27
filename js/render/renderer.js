// js/render/renderer.js
// (facade, split into modules in Wave R1)

import { render_image_tem, draw_atoms, draw_bonds } from './render_tem.js';
import { render_afm_like } from './render_afm.js';
import { render_stm_like } from './render_stm.js';
import { compute_scaled_coordinates } from './render_common.js';

export { compute_scaled_coordinates, draw_atoms, draw_bonds };

export function render_image(atoms, opts = {}) {
    const m = String((opts && opts.mode) ? opts.mode : 'TEM').trim().toUpperCase();
    if (m === 'AFM') return render_afm_like(atoms, opts);
    if (m === 'STM') return render_stm_like(atoms, opts);
    return render_image_tem(atoms, opts);
}
