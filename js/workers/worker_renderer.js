// js/worker_renderer.js
import { render_image } from '../render/renderer.js';

self.onmessage = (ev) => {
    const { atoms, opts } = ev.data;
    const out = render_image(atoms, opts); // Uint8ClampedArray
    self.postMessage({ frame: out, w: opts.img_size[1], h: opts.img_size[0] }, [out.buffer]);
};
