// js/app/main.js
// UI + state + rebuild/render pipeline (simplified)
//
// IMPORTANT:
//  - main.js MUST NOT construct ImageData from a grayscale buffer.
//    renderer.js already converts grayscale -> RGBA and draws into canvas.
//
// UI (current):
//  - CIF/CSV are loaded ONLY via file input / drag&drop.
//  - No "liquid/solid" modes.

import { build_system_from_input } from './phases.js';
import { render_image } from '../render/renderer.js';
import { make_view_state } from './camera.js';
import { GifRecorder } from '../export/gif_recorder.js';

export async function interactive_em_image(smiles_text, _unused) {
    function el(id) { return document.getElementById(id); }

    // ---- DOM ----
    var cvs = el('canvas');
    if (!cvs) throw new Error('Canvas #canvas not found');
    var ctx = cvs.getContext('2d', { willReadFrequently: true });

    var lblTitle = el('lbl-title');

    // Controls
    var rngZoom = el('rng-zoom');
    var rngContrast = el('rng-contrast');
    var rngBlur = el('rng-blur');
    var rngBg = el('rng-bg');
    var lblBg = el('lbl-bg');
    var rngNoise = el('rng-noise');
    var rngFocus = el('rng-focus');
    var rngDof = el('rng-dof');
    var rngBwidth = el('rng-bwidth');
    var rngBamp = el('rng-bamp');

    var cbInvert = el('cb-invert');
    var cbNoise = el('cb-noise');
    var cbBonds = el('cb-bonds');
    var cbHide = el('cb-hidefront');
    var cbScale = el('cb-scale');

    // Bonds label (AFM UI rename only)
    var lblBonds = el('lblBonds');
    var bondsLabelDefault = (lblBonds && typeof lblBonds.textContent === 'string') ? lblBonds.textContent : '';

    // Microscopy mode presets (UI-only in Wave 1)
    var selMode = el('sel-mode');
    var lblModeDesc = el('lbl-mode-desc');
    var cbScanlines = el('cb-scanlines');
    var rowScanlines = el('row-scanlines');
    var rngTip = el('rng-tip');
    var rowTip = el('row-tip');

    var tbTiles = el('tb-tiles');
    var tbSmiles = el('tb-smiles');
    var tbW = el('tb-w');
    var tbH = el('tb-h');

    var fileInput = el('file-input');

    // Buttons
    var btnExport = el('btn-export');

    // GIF buttons
    var btnGifStart = el('btn-gif-start');
    var btnGifStop = el('btn-gif-stop');
    var btnGifDownload = el('btn-gif-download');
    var lblGifStatus = el('lbl-gif-status');

    // View (stage) buttons
    var btnViewUp = el('btn-view-up');
    var btnViewDown = el('btn-view-down');
    var btnViewLeft = el('btn-view-left');
    var btnViewRight = el('btn-view-right');
    var btnViewRotL = el('btn-view-rotl');
    var btnViewRotR = el('btn-view-rotr');
    var btnViewReset = el('btn-view-reset');
    var btnModelX90 = el('btn-model-x90');
    var btnModelY90 = el('btn-model-y90');


    // ---- State ----
    // View state (UI-level; used by provider selection and renderer projection)
    // pan_px: [dx, dy] in pixels; rot_deg: rotation around Z in degrees (XY plane)
    var view = {
        pan_px: [0, 0],
        rot_deg: 0,
        tilt_x_deg: 0,
        tilt_y_deg: 0,
        angstroms_per_pixel: 0.1,
        img_size: [400, 400],
        center_A: [0, 0, 0],   // anchor in Å (XY; Z kept for completeness)
        canvas_w: 400,
        canvas_h: 400
    };

    const PAN_STEP_PX = 20;
    const ROT_STEP_DEG = 5;
    const TILT_STEP_DEG = 90;

    // UI microscopy mode preset (does NOT change renderer math in Wave 1)
    var mode = 'TEM';

    // Provider-based state (lazy-ready)
    var provider = null;
    var providerMeta = null;
    var title = '';

    // last chosen source
    var active_source = 'smiles'; // 'smiles' | 'file'

    // currently loaded file (if any)
    var file_state = {
        kind: null,  // 'cif' | 'csv' | 'xyz' | 'json' | 'poscar' | null
        name: '',
        url: null,
        isBlob: false
    };


    // ---- GIF recording ----
    // Recording is session-based: Start begins (or resumes) a session, Stop pauses it,
    // Download finalizes and unlocks size controls.
    var gifRec = new GifRecorder(cvs, ctx, { maxColors: 256, repeat: 0, maxFrames: 900 });

    // Lock W/H while GIF session is active (even paused). Unlocked only after download.
    var sizeLock = { active: false, w: 0, h: 0 };

    function setSizeLock(on) {
        sizeLock.active = !!on;
        if (sizeLock.active) {
            sizeLock.w = cvs.width | 0;
            sizeLock.h = cvs.height | 0;
        }
        if (tbW) tbW.disabled = sizeLock.active;
        if (tbH) tbH.disabled = sizeLock.active;
    }

    function updateGifUI() {
        if (!btnGifStart || !btnGifStop || !btnGifDownload) return;

        var session = gifRec.isSessionActive();
        var recording = gifRec.isRecording();

        btnGifStart.disabled = false;
        btnGifStop.disabled = !session;
        btnGifDownload.disabled = !session;

        if (lblGifStatus) {
            if (!session) {
                lblGifStatus.textContent = '';
            } else {
                var n = gifRec.getFrameCount();
                lblGifStatus.textContent = recording
                    ? (tr('gifRecPrefix', 'GIF: recording… frames=') + n)
                    : (tr('gifPausePrefix', 'GIF: paused. frames=') + n + tr('gifPauseSuffix', ' (Download to finish)'));
}
        }
    }

    
    // i18n helper (uses globals set by translations loader in index.html/help.html)
    function tr(key, fallback) {
        try {
            var d = window.__I18N_DICT;
            if (d && typeof d[key] === 'string') return d[key];
        } catch (e) {}
        return fallback;
    }

// ---- utils ----
    function trimStr(s) { return (s || '').trim(); }

    function clampInt(n, lo, hi, fallback) {
        var v = parseInt(n, 10);
        if (!Number.isFinite(v)) v = fallback;
        v = Math.min(hi, Math.max(lo, v | 0));
        return v;
    }

    function getBackgroundGray() {
        var bg = rngBg ? parseInt(rngBg.value, 10) : 127;
        if (!Number.isFinite(bg)) bg = 127;
        bg = Math.min(255, Math.max(0, bg | 0));
        if (lblBg) lblBg.textContent = String(bg);
        return bg;
    }

    // ---- microscopy mode presets (Wave 1: UI/UX only) ----
    function normMode(m) {
        m = String(m || '').trim().toUpperCase();
        if (m === 'TEM' || m === 'STM' || m === 'AFM') return m;
        return 'TEM';
    }

    function clampRangeValue(rng, v) {
        if (!rng) return null;
        var lo = parseFloat(rng.min);
        var hi = parseFloat(rng.max);
        var x = parseFloat(v);
        if (!Number.isFinite(x)) x = parseFloat(rng.value);
        if (Number.isFinite(lo)) x = Math.max(lo, x);
        if (Number.isFinite(hi)) x = Math.min(hi, x);
        // keep string formatting (avoid scientific notation in UI)
        return String(x);
    }

    function setHidden(node, hidden) {
        if (!node) return;
        try { node.classList.toggle('ui_hidden', !!hidden); } catch (e) {
            node.style.display = hidden ? 'none' : '';
        }
    }

    function setWrapDisabled(input, disabled) {
        if (!input) return;
        input.disabled = !!disabled;

        var wrap = null;
        try {
            wrap = input.closest('label') || input.closest('.range_row');
        } catch (e) { wrap = null; }

        if (wrap) {
            try { wrap.classList.toggle('ui_disabled', !!disabled); } catch (e) {}
        } else {
            // fallback: direct style
            try { input.style.opacity = disabled ? '0.5' : ''; } catch (e) {}
        }
    }

    function updateModeDescText() {
        if (!lblModeDesc) return;
        if (mode === 'STM') lblModeDesc.textContent = tr('modeDescSTM', 'STM — поверхневий режим, DoF вимкнено.');
        else if (mode === 'AFM') lblModeDesc.textContent = tr('modeDescAFM', 'AFM — поверхневий режим, DoF вимкнено.');
        else lblModeDesc.textContent = tr('modeDescTEM', 'TEM — “проєкційний контраст”, STM/AFM — “поверхневі режими, DoF вимкнено”.');
    }

    function applyModePreset(newMode, cfg) {
        cfg = cfg || {};
        var save = (cfg.save !== false);
        var rerender = (cfg.rerender !== false);

        mode = normMode(newMode);
        if (selMode) selMode.value = mode;

        if (save) {
            try { localStorage.setItem('scigentem_mode', mode); } catch (e) {}
        }

        // Visibility for mode-specific controls
        setHidden(rowScanlines, mode !== 'STM');
        setHidden(rowTip, mode !== 'AFM');

        // Global policies (Wave 1.1):
        //  - DoF is always reset to 0 on mode switch, but NEVER disabled.
        //  - Focus Z + Hide front are NEVER disabled.
        if (rngDof) rngDof.value = clampRangeValue(rngDof, 0.0);
        setWrapDisabled(rngDof, false);
        setWrapDisabled(cbHide, false);

        // AFM only: rename bonds label in UI (no semantic change)
        if (lblBonds) {
            if (mode === 'AFM') lblBonds.textContent = 'Bond contrast (AFM)';
            else lblBonds.textContent = bondsLabelDefault || 'Зв\'язки';
        }

        // TEM: bonds OFF + disabled, DoF reset to 0 (but allowed), blur/noise moderate
        if (mode === 'TEM') {
            if (cbBonds) cbBonds.checked = false;
            setWrapDisabled(cbBonds, true);

            if (cbInvert) cbInvert.checked = false;
            if (rngBlur) rngBlur.value = clampRangeValue(rngBlur, 1.0);

            if (cbNoise) cbNoise.checked = true;
            if (rngNoise) rngNoise.value = clampRangeValue(rngNoise, 3.0);

            if (rngContrast) rngContrast.value = clampRangeValue(rngContrast, 1.35);

            if (cbScanlines) cbScanlines.checked = false;
        }

        // STM: bonds OFF + disabled, DoF reset to 0 (but allowed), higher contrast, smaller blur
        if (mode === 'STM') {
            if (cbBonds) cbBonds.checked = false;
            setWrapDisabled(cbBonds, true);

            if (cbInvert) cbInvert.checked = false;
            if (rngBlur) rngBlur.value = clampRangeValue(rngBlur, 0.4);

            if (cbNoise) cbNoise.checked = true;
            if (rngNoise) rngNoise.value = clampRangeValue(rngNoise, 1.6);

            if (rngContrast) rngContrast.value = clampRangeValue(rngContrast, 2.2);

            if (cbScanlines) cbScanlines.checked = true;
            setWrapDisabled(cbScanlines, false);
            setWrapDisabled(rngTip, true);
        }

        // AFM: bonds ON by default, DoF reset to 0 (but allowed), moderate blur/noise, mid-high contrast
        if (mode === 'AFM') {
            if (cbBonds) cbBonds.checked = true;
            setWrapDisabled(cbBonds, false);

            if (cbInvert) cbInvert.checked = false;
            if (rngBlur) rngBlur.value = clampRangeValue(rngBlur, 0.8);

            if (cbNoise) cbNoise.checked = true;
            if (rngNoise) rngNoise.value = clampRangeValue(rngNoise, 1.2);

            if (rngContrast) rngContrast.value = clampRangeValue(rngContrast, 1.7);

            if (rngTip) rngTip.value = clampRangeValue(rngTip, 0.5);
            setWrapDisabled(rngTip, false);
            setWrapDisabled(cbScanlines, true);
        }

        // Keep background label in sync
        if (rngBg) getBackgroundGray();

        updateModeDescText();
        if (rerender) scheduleRender();
    }



    function clearCanvas(bg) {
        ctx.save();
        ctx.fillStyle = 'rgb(' + bg + ',' + bg + ',' + bg + ')';
        ctx.fillRect(0, 0, cvs.width, cvs.height);
        ctx.restore();
    }

    function setTitle(t) {
        title = t || '';
        if (lblTitle) lblTitle.textContent = 'Simulated EM — ' + title;
    }

    function fileKindFromName(name) {
        var low = String(name || '').toLowerCase();
        // tolerate things like "file.cif;..." or "file.cif?x" in the #fragment
        if (low.indexOf('.cif') >= 0) return 'cif';
        if (low.indexOf('.csv') >= 0) return 'csv';
        if (low.indexOf('.xyz') >= 0) return 'xyz';
        if (low.indexOf('.json') >= 0) return 'json';
        if (low.indexOf('.mol') >= 0) return 'mol';
        if (low.indexOf('.poscar') >= 0) return 'poscar';
        if (low.indexOf('.contcar') >= 0) return 'poscar';
        if (low.indexOf('.vasp') >= 0) return 'poscar';
        // common VASP names without extension
        if (low === 'poscar' || low === 'contcar') return 'poscar';
        return null;
    }

    function revokeFileUrl() {
        if (file_state.isBlob && file_state.url) {
            try { URL.revokeObjectURL(file_state.url); } catch (e) { }
        }
        file_state.kind = null;
        file_state.name = '';
        file_state.url = null;
        file_state.isBlob = false;
    }

    function parseTilesSize() {
        if (!tbTiles) return null;
        var t = trimStr(tbTiles.value);
        if (!t) return null;
        if (t.toLowerCase() === 'auto') return 'auto';
        return t.replace(/×/g, 'x');
    }

    function applyPeriodicSizeToSpec(spec, sizeStr) {
        if (!sizeStr) return spec;

        var s = String(spec || '').trim();
        if (!s) return s;

        var semi = s.indexOf(';');
        if (semi < 0) return s + ';size=' + sizeStr;

        var base = s.slice(0, semi).trim();
        var rest = s.slice(semi + 1).trim();
        var items = rest ? rest.split(';') : [];

        var kept = [];
        for (var i = 0; i < items.length; i++) {
            var it = trimStr(items[i]);
            if (!it) continue;
            var itLow = it.toLowerCase();
            if (itLow.indexOf('size=') === 0) continue;
            kept.push(it);
        }
        kept.push('size=' + sizeStr);

        return base + ';' + kept.join(';');
    }

    function isPeriodicFileKind(kind, name) {
        if (kind === 'cif' || kind === 'json' || kind === 'poscar') return true;
        var low = String(name || '').toLowerCase();
        if (low === 'poscar' || low === 'contcar') return true;
        return false;
    }

    // ---- build + render ----
    async function buildCurrentSystem() {
        var spec = '';
        var kind = null;

        if (active_source === 'file' && file_state.url) {
            // extension is encoded in #fragment so phases can detect .cif/.csv
            spec = file_state.url + '#' + file_state.name;
            kind = file_state.kind;

            if (isPeriodicFileKind(kind, file_state.name)) {
                var tiles = parseTilesSize();
                if (tiles) spec = applyPeriodicSizeToSpec(spec, tiles);
            }
        } else {
            active_source = 'smiles';
            spec = (tbSmiles ? trimStr(tbSmiles.value) : '') || 'O';
            kind = 'smiles';
        }

        provider = await build_system_from_input(spec);
        providerMeta = (provider && typeof provider.getMeta === 'function') ? provider.getMeta() : null;
        if (providerMeta && Array.isArray(providerMeta.center_A)) {
            // Keep a stable anchor for pan/rot and for future lazy-tiling selection.
            view.center_A = providerMeta.center_A;
        }
        var uiTitle = (active_source === 'file' && file_state && file_state.name) ? file_state.name
            : ((providerMeta && providerMeta.title) ? providerMeta.title : (kind || ''));
        setTitle(uiTitle);
    }

    async function renderCurrent() {
        var w = clampInt(tbW ? tbW.value : 400, 128, 4096, 400);
        var h = clampInt(tbH ? tbH.value : 400, 128, 4096, 400);
        if (tbW) tbW.value = String(w);
        if (tbH) tbH.value = String(h);

        if (cvs.width !== w) cvs.width = w;
        if (cvs.height !== h) cvs.height = h;

        if (!provider || !providerMeta) {
            clearCanvas(getBackgroundGray());
            return;
        }

        // keep viewState in sync with current render params
        view.canvas_w = w;
        view.canvas_h = h;
        view.img_size = [h, w];
        view.angstroms_per_pixel = rngZoom ? parseFloat(rngZoom.value) : view.angstroms_per_pixel;

        // build viewState early (camera affects BOTH selection and projection)
        var viewState = make_view_state({
            img_size: view.img_size,
            angstroms_per_pixel: view.angstroms_per_pixel,
            pan_px: view.pan_px,
            rotZ_rad: (view.rot_deg * Math.PI) / 180.0,
            center_A: (view.center_A || [0, 0, 0]),
            center_mode: 'bbox'
        });
        // optional compatibility fields (provider/camera helpers may use them later)
        viewState.canvas_w = view.canvas_w;
        viewState.canvas_h = view.canvas_h;
        viewState.rot_deg = view.rot_deg;

        // model tilts (rotate the structure itself in 3D before projection)
        viewState.rotX_rad = (view.tilt_x_deg * Math.PI) / 180.0;
        viewState.rotY_rad = (view.tilt_y_deg * Math.PI) / 180.0;
        viewState.tilt_x_deg = view.tilt_x_deg;
        viewState.tilt_y_deg = view.tilt_y_deg;

        var provView = provider.getView(viewState, { needBonds: (cbBonds ? !!cbBonds.checked : true) });
        var atoms = (provView && provView.atomsView) ? provView.atomsView : [];
        var bonds = (provView && ('bondsView' in provView)) ? provView.bondsView : null;
        var usedCamera = (provView && provView.usedCamera) ? provView.usedCamera : viewState;
        if (usedCamera) {
            usedCamera.canvas_w = view.canvas_w;
            usedCamera.canvas_h = view.canvas_h;
            usedCamera.rot_deg = view.rot_deg;
            usedCamera.rotX_rad = viewState.rotX_rad;
            usedCamera.rotY_rad = viewState.rotY_rad;
            usedCamera.tilt_x_deg = view.tilt_x_deg;
            usedCamera.tilt_y_deg = view.tilt_y_deg;
        }

        if (!atoms || atoms.length === 0) {
            clearCanvas(getBackgroundGray());
            return;
        }

        // focal_z from atoms
        var focus_norm = rngFocus ? parseFloat(rngFocus.value) : 0.5;
        if (!Number.isFinite(focus_norm)) focus_norm = 0.5;

        // focal_z from atoms (respect model tilt X/Y if used)
        var cxA = (usedCamera && usedCamera.center_A) ? (usedCamera.center_A[0] || 0) : 0;
        var cyA = (usedCamera && usedCamera.center_A) ? (usedCamera.center_A[1] || 0) : 0;
        var czA = (usedCamera && usedCamera.center_A) ? (usedCamera.center_A[2] || 0) : 0;
        var ax = (usedCamera && Number.isFinite(usedCamera.rotX_rad)) ? usedCamera.rotX_rad : 0;
        var ay = (usedCamera && Number.isFinite(usedCamera.rotY_rad)) ? usedCamera.rotY_rad : 0;

        function zAfterModelTilt(a) {
            var x = (a.x || 0) - cxA;
            var y = (a.y || 0) - cyA;
            var z = ((a.z != null ? a.z : 0) - czA);

            if (ax) {
                var c = Math.cos(ax), s = Math.sin(ax);
                var y1 = y * c - z * s;
                var z1 = y * s + z * c;
                y = y1;
                z = z1;
            }
            if (ay) {
                var c2 = Math.cos(ay), s2 = Math.sin(ay);
                var x2 = x * c2 + z * s2;
                var z2 = -x * s2 + z * c2;
                x = x2;
                z = z2;
            }
            return z + czA;
        }

        var useTilt = (Math.abs(ax) > 1e-12) || (Math.abs(ay) > 1e-12);
        var zmin = useTilt ? zAfterModelTilt(atoms[0]) : ((atoms[0].z != null) ? atoms[0].z : 0);
        var zmax = zmin;
        for (var i = 1; i < atoms.length; i++) {
            var zc = useTilt ? zAfterModelTilt(atoms[i]) : ((atoms[i].z != null) ? atoms[i].z : 0);
            if (zc < zmin) zmin = zc;
            if (zc > zmax) zmax = zc;
        }
        var focal_z = zmin + focus_norm * (zmax - zmin);

        var opts = {
            bonds: bonds,
            img_size: [h, w],
            angstroms_per_pixel: viewState.angstroms_per_pixel,
            blur_sigma: rngBlur ? parseFloat(rngBlur.value) : 1.0,
            background_gray: getBackgroundGray(),
            invert: cbInvert ? !!cbInvert.checked : false,
            noise_stddev: (cbNoise && cbNoise.checked && rngNoise) ? parseFloat(rngNoise.value) : 0.0,
            contrast: rngContrast ? parseFloat(rngContrast.value) : 1.0,
            compose_mode: 'sum',
            draw_bonds_flag: cbBonds ? !!cbBonds.checked : true,
            camera: usedCamera,
            bond_wave_width_px: rngBwidth ? parseFloat(rngBwidth.value) : 6,
            bond_wave_amplitude: rngBamp ? parseFloat(rngBamp.value) : 0.4,
            low_clip: (el('tb-clip-lo') && el('tb-clip-lo').value !== '') ? parseFloat(el('tb-clip-lo').value) : null,
            high_clip: (el('tb-clip-hi') && el('tb-clip-hi').value !== '') ? parseFloat(el('tb-clip-hi').value) : null,
            focal_z: focal_z,
            dof_strength: rngDof ? parseFloat(rngDof.value) : 0.0,
            hide_front: cbHide ? !!cbHide.checked : false,
            show_scale_bar: cbScale ? !!cbScale.checked : false,
            scale_bar_corner: 'bl',
            scale_bar_margin_px: 12,
            canvasCtx: ctx,

            // UI-only (Wave 1): renderer may ignore these for now
            mode: mode,
            scanlines: (cbScanlines ? !!cbScanlines.checked : false),
            tip_sharpness: (rngTip ? parseFloat(rngTip.value) : 0.5)
        };

        // renderer.js draws into canvasCtx
        render_image(atoms, opts);

        // If recording, capture the current canvas as a GIF frame
        if (gifRec && gifRec.isRecording && gifRec.isRecording()) {
            gifRec.captureFrame();
        }
        updateGifUI();
    }

    var rafPending = false;
    function scheduleRender() {
        if (rafPending) return;
        rafPending = true;
        requestAnimationFrame(function () {
            rafPending = false;
            renderCurrent();
        });
    }

    var rebuildPending = false;
    var rebuildQueued = false;
    async function rebuildAndRender() {
        // Coalesce rebuild requests instead of dropping them.
        // This prevents the "typed SMILES -> temporary invalid -> fallback O -> never updates" issue
        // when RDKit build is still running.
        if (rebuildPending) {
            rebuildQueued = true;
            return;
        }

        rebuildPending = true;
        try {
            await buildCurrentSystem();
        } catch (e) {
            // Keep UI responsive even on parse errors.
            // (Minimal error handling; no debug spam.)
            provider = null;
            providerMeta = null;
            setTitle('ERROR');
        } finally {
            rebuildPending = false;
        }

        await renderCurrent();

        if (rebuildQueued) {
            rebuildQueued = false;
            setTimeout(function () { rebuildAndRender(); }, 0);
        }
    }

    // ---- file load ----
    async function handleLocalFile(file) {
        if (!file) return;
        var kind = fileKindFromName(file.name);
        // Allow extensionless POSCAR/CONTCAR and content-sniffing in phases.js (blob urls).
        // Keep a light safety net: ignore huge binary files.
        if (!kind && (file.size | 0) > (20 * 1024 * 1024)) return; // 20 MB

        revokeFileUrl();

        file_state.kind = kind;
        file_state.name = file.name;
        file_state.url = URL.createObjectURL(file);
        file_state.isBlob = true;

        active_source = 'file';
        await rebuildAndRender();
    }

    // ---- events ----

    // Render-only controls
    [rngZoom, rngContrast, rngBlur, rngBg, rngNoise, rngFocus, rngDof, rngBwidth, rngBamp, rngTip].forEach(function (x) {
        if (!x) return;
        x.oninput = function () {
            if (x === rngBg) getBackgroundGray();
            scheduleRender();
        };
    });

    [cbInvert, cbNoise, cbBonds, cbHide, cbScale, cbScanlines].forEach(function (x) {
        if (x) x.onchange = scheduleRender;
    });

    if (selMode) {
        selMode.onchange = function () {
            applyModePreset(selMode.value, { save: true, rerender: true });
        };
    }

    if (tbW) tbW.onchange = function () {
        if (sizeLock && sizeLock.active) { tbW.value = String(sizeLock.w); return; }
        scheduleRender();
    };
    if (tbH) tbH.onchange = function () {
        if (sizeLock && sizeLock.active) { tbH.value = String(sizeLock.h); return; }
        scheduleRender();
    };

    if (tbTiles) {
        tbTiles.onchange = function () {
            if (active_source === 'file' && isPeriodicFileKind(file_state.kind, file_state.name)) rebuildAndRender();
        };
    }

    if (tbSmiles) {
        tbSmiles.value = trimStr(smiles_text);

        // Switching from file -> SMILES must be immediate.
        // We rebuild after a short debounce to avoid running RDKit on every keystroke.
        var smilesTimer = null;
        function scheduleSmilesRebuild() {
            active_source = 'smiles';
            if (smilesTimer) clearTimeout(smilesTimer);
            smilesTimer = setTimeout(function () {
                smilesTimer = null;
                rebuildAndRender();
            }, 250);
        }

        tbSmiles.oninput = scheduleSmilesRebuild;
        tbSmiles.onchange = function () {
            active_source = 'smiles';
            rebuildAndRender();
        };

        // Make plain Enter apply.
        // - For <textarea>: Enter applies; Shift+Enter inserts newline.
        // - For <input>: Enter applies.
        tbSmiles.addEventListener('keydown', function (ev) {
            if (ev.key !== 'Enter') return;

            var isTextarea = String(tbSmiles.tagName || '').toUpperCase() === 'TEXTAREA';
            if (isTextarea && ev.shiftKey) return; // allow newline

            ev.preventDefault();
            active_source = 'smiles';
            rebuildAndRender();
        });
    }

    if (fileInput) {    
        fileInput.onchange = function () {
            var files = fileInput.files;
            if (files && files.length) handleLocalFile(files[0]);
            fileInput.value = '';
        };
    }

    // Global drag&drop (drop anywhere)
    document.addEventListener('dragover', function (ev) { ev.preventDefault(); });
    document.addEventListener('drop', function (ev) {
        ev.preventDefault();
        var dt = ev.dataTransfer;
        if (!dt || !dt.files || !dt.files.length) return;
        handleLocalFile(dt.files[0]);
    });

    if (btnExport) {
        btnExport.onclick = function () {
            var a = document.createElement('a');
            a.download = 'tem.png';
            a.href = cvs.toDataURL('image/png');
            a.click();
        };
    }


    // GIF controls
    if (btnGifStart) {
        btnGifStart.onclick = async function () {
            // Start/resume session and lock size immediately
            if (gifRec && !gifRec.isSessionActive()) {
                // ensure current render size is applied before locking
                await renderCurrent();
                setSizeLock(true);
            }
            if (gifRec) await gifRec.startOrResume();
            updateGifUI();
            // force a render so the first frame is captured
            scheduleRender();
        };
    }

    if (btnGifStop) {
        btnGifStop.onclick = function () {
            if (gifRec) gifRec.pause();
            updateGifUI();
        };
    }

    if (btnGifDownload) {
        btnGifDownload.onclick = function () {
            if (!gifRec) return;
            var blob = gifRec.finishToBlob();
            if (!blob) {
                // no frames; unlock and reset
                setSizeLock(false);
                updateGifUI();
                return;
            }

            var a = document.createElement('a');
            a.download = 'tem.gif';
            a.href = URL.createObjectURL(blob);
            a.click();

            setTimeout(function () {
                try { URL.revokeObjectURL(a.href); } catch (e) {}
            }, 1000);

            // unlock size only after download finalizes
            setSizeLock(false);
            updateGifUI();
        };
    }

    // View controls: pan (px) + rotation (deg) — rerender only (NO rebuild/parsing)
    function updateViewAndRender(fn) {
        if (typeof fn === 'function') fn();
        scheduleRender();
    }

    if (btnViewUp) btnViewUp.onclick = function () { updateViewAndRender(function () { view.pan_px[1] -= PAN_STEP_PX; }); };
    if (btnViewDown) btnViewDown.onclick = function () { updateViewAndRender(function () { view.pan_px[1] += PAN_STEP_PX; }); };
    if (btnViewLeft) btnViewLeft.onclick = function () { updateViewAndRender(function () { view.pan_px[0] -= PAN_STEP_PX; }); };
    if (btnViewRight) btnViewRight.onclick = function () { updateViewAndRender(function () { view.pan_px[0] += PAN_STEP_PX; }); };
    if (btnViewRotL) btnViewRotL.onclick = function () { updateViewAndRender(function () { view.rot_deg -= ROT_STEP_DEG; }); };
    if (btnViewRotR) btnViewRotR.onclick = function () { updateViewAndRender(function () { view.rot_deg += ROT_STEP_DEG; }); };
    if (btnViewReset) btnViewReset.onclick = function () { updateViewAndRender(function () { view.pan_px = [0, 0]; view.rot_deg = 0; }); };

    // Model tilt: rotate structure itself by +90° around X or Y
    if (btnModelX90) btnModelX90.onclick = function () {
        updateViewAndRender(function () {
            view.tilt_x_deg = (view.tilt_x_deg + TILT_STEP_DEG) % 360;
        });
    };
    if (btnModelY90) btnModelY90.onclick = function () {
        updateViewAndRender(function () {
            view.tilt_y_deg = (view.tilt_y_deg + TILT_STEP_DEG) % 360;
        });
    };

    // refresh dynamic labels on language switch
    try {
        document.addEventListener('i18n-updated', function () { updateGifUI(); updateModeDescText(); });
    } catch (e) {}

    // ---- init ----
    // Restore microscopy mode (presets only affect UI + render params, not parsing/build)
    try {
        var storedMode = normMode(localStorage.getItem('scigentem_mode'));
        applyModePreset(storedMode, { save: false, rerender: false });
    } catch (e) {
        applyModePreset('TEM', { save: false, rerender: false });
    }

    active_source = 'smiles';
    getBackgroundGray();
    updateGifUI();
    await rebuildAndRender();
}
