// js/export/gif_recorder.js
// Client-side GIF recorder for Canvas 2D using gifenc (ESM via CDN).
// Captures frames from canvas ctx.getImageData() and encodes to GIF bytes.
//
// Notes:
// - The recorder captures frames when captureFrame() is called (typically after each render).
// - Size is fixed for a session; changing canvas size mid-session is not supported.

const GIFENC_URL = 'https://unpkg.com/gifenc@1.0.3?module';

let _gifencPromise = null;
async function loadGifenc() {
  if (!_gifencPromise) _gifencPromise = import(GIFENC_URL);
  return _gifencPromise;
}

export class GifRecorder {
  constructor(canvas, ctx, opts = {}) {
    this.canvas = canvas;
    this.ctx = ctx;

    this.maxColors = opts.maxColors ?? 256;
    this.repeat = opts.repeat ?? 0;     // 0=forever, -1=once
    this.maxFrames = opts.maxFrames ?? 900;

    this.sessionActive = false; // started at least once, not finalized yet
    this.recording = false;     // capturing enabled
    this.frameCount = 0;

    this._gif = null;
    this._quantize = null;
    this._applyPalette = null;

    this._lastTs = 0;
    this._status = 'idle';
  }

  isSessionActive() { return this.sessionActive; }
  isRecording() { return this.sessionActive && this.recording; }

  getStatus() { return this._status; }
  getFrameCount() { return this.frameCount; }

  async startOrResume() {
    const { GIFEncoder, quantize, applyPalette } = await loadGifenc();

    if (!this.sessionActive) {
      this._gif = GIFEncoder();
      this._quantize = quantize;
      this._applyPalette = applyPalette;

      this.sessionActive = true;
      this.frameCount = 0;
      this._lastTs = performance.now();
    }

    this.recording = true;
    this._status = 'recording';
  }

  pause() {
    if (!this.sessionActive) return;
    this.recording = false;
    this._status = 'paused';
  }

  reset() {
    this.sessionActive = false;
    this.recording = false;
    this.frameCount = 0;
    this._lastTs = 0;
    this._gif = null;
    this._quantize = null;
    this._applyPalette = null;
    this._status = 'idle';
  }

  // Call right AFTER rendering into the canvas.
  captureFrame() {
    if (!this.isRecording()) return { ok: false, reason: 'not-recording' };
    if (!this._gif || !this._quantize || !this._applyPalette) return { ok: false, reason: 'not-ready' };
    if (this.frameCount >= this.maxFrames) {
      this.recording = false;
      this._status = 'paused(maxFrames)';
      return { ok: false, reason: 'maxFrames' };
    }

    const w = this.canvas.width | 0;
    const h = this.canvas.height | 0;
    if (w <= 0 || h <= 0) return { ok: false, reason: 'bad-size' };

    let img;
    try {
      img = this.ctx.getImageData(0, 0, w, h);
    } catch (e) {
      this.recording = false;
      this._status = 'error(getImageData)';
      return { ok: false, reason: 'getImageData' };
    }

    const now = performance.now();
    let delay = now - (this._lastTs || now);
    this._lastTs = now;

    // clamp to sane range
    if (!Number.isFinite(delay)) delay = 100;
    delay = Math.max(20, Math.min(2000, delay));

    const rgba = img.data; // Uint8ClampedArray RGBA
    const palette = this._quantize(rgba, this.maxColors);
    const index = this._applyPalette(rgba, palette);

    this._gif.writeFrame(index, w, h, {
      palette,
      delay: Math.round(delay),
      repeat: (this.frameCount === 0) ? this.repeat : undefined
    });

    this.frameCount++;
    return { ok: true };
  }

  finishToBlob() {
    if (!this.sessionActive || !this._gif) return null;
    if (this.frameCount <= 0) {
      this.reset();
      return null;
    }

    this.recording = false;
    this._gif.finish();
    const bytes = this._gif.bytes();
    const blob = new Blob([bytes], { type: 'image/gif' });

    // reset session
    this.reset();
    this._status = 'finished';

    return blob;
  }
}
