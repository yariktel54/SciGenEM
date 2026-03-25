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

import { build_system_from_input } from "./phases.js";
import { render_image } from "../render/renderer.js";
import {
  make_view_state,
  panBy,
  resetPan,
  rotateBy,
  resetRotation,
} from "./camera.js";
import { GifRecorder } from "../export/gif_recorder.js";

export async function interactive_em_image(smiles_text, _unused) {
  function el(id) {
    return document.getElementById(id);
  }

  // ---- DOM ----
  var cvs = el("canvas");
  if (!cvs) throw new Error("Canvas #canvas not found");
  var ctx = cvs.getContext("2d", { willReadFrequently: true });

  var lblTitle = el("lbl-title");

  // Controls
  var rngZoom = el("rng-zoom");
  var rngContrast = el("rng-contrast");
  var rngBlur = el("rng-blur");
  var rngBg = el("rng-bg");
  var lblBg = el("lbl-bg");
  var rngNoise = el("rng-noise");
  var rngFocus = el("rng-focus");
  var rngDof = el("rng-dof");
  var rngBwidth = el("rng-bwidth");
  var rngBamp = el("rng-bamp");
  // Brightness clip (min/max) boxes (now under "Additional settings")
  var tbClipLo = el("tb-clip-lo");
  var tbClipHi = el("tb-clip-hi");

  var cbInvert = el("cb-invert");
  var cbNoise = el("cb-noise");
  var cbBonds = el("cb-bonds");
  var cbHide = el("cb-hidefront");
  var cbScale = el("cb-scale");
  var cbTipHighlight = el("cb-tip-highlight");
  var cbTipProject = el("cb-tip-project");

  // Bonds label (UI text; dynamically changes per mode and language)
  var lblBonds = el("lblBonds");

  // Microscopy mode presets (UI-only in Wave 1)
  var selMode = el("sel-mode");
  var lblModeDesc = el("lbl-mode-desc");
  var cbScanlines = el("cb-scanlines");
  var rowScanlines = el("row-scanlines");
  var rngTip = el("rng-tip");
  var rowTip = el("row-tip");

  var tbNx = el("tb-nx");
  var tbNy = el("tb-ny");
  var tbNz = el("tb-nz");
  var tbSmiles = el("tb-smiles");
  var lblSmilesErr = el("smiles-error");
  var selSamples = el("sel-samples");
  var lblSamplesStatus = el("lblSamplesStatus");
  var btnResetSettings = el("resetSettings");
  var tbW = el("tb-w");
  var tbH = el("tb-h");

  var fileInput = el("file-input");
  var sceneList = el("scene-list");
  var sceneEmpty = el("scene-empty");
  var sceneStats = el("scene-stats");

  var lblMoveValues = el("lbl-move-values");
  var lblRotateValues = el("lbl-rotate-values");
  var btnMoveUp = el("btn-pan-up");
  var btnMoveDown = el("btn-pan-down");
  var btnMoveLeft = el("btn-pan-left");
  var btnMoveRight = el("btn-pan-right");
  var btnMoveHome = el("btn-pan-home");
  var btnMoveZp = el("btn-move-zp");
  var btnMoveZm = el("btn-move-zm");
  var btnRotateXp = el("btn-rot-xp");
  var btnRotateXm = el("btn-rot-xm");
  var btnRotateYp = el("btn-rot-yp");
  var btnRotateYm = el("btn-rot-ym");
  var btnRotateHome = el("btn-rot-home");
  var btnRotateZp = el("btn-rot-zp");
  var btnRotateZm = el("btn-rot-zm");
  var btnResetObjectEdits = el("btn-reset-object-edits");
  var btnResetAllObjectEdits = el("btn-reset-all-object-edits");
  var clipZDetails = el("object-advanced-details");
  var cbClipZEnabled = el("cb-clipz-enabled");
  var rngClipZMin = el("rng-clipz-min");
  var rngClipZMax = el("rng-clipz-max");
  var lblClipZValues = el("lbl-clipz-values");
  var btnClipZReset = el("btn-clipz-reset");

  // Element overrides (per-element multipliers)
  var elemOvTableWrap = el("elem-override-table");
  var elemOvDetails = el("elem-override-details");
  var btnElemOvReset = el("btn-elem-ov-reset");
  var btnElemOvExport = el("btn-elem-ov-export");
  var btnElemOvImport = el("btn-elem-ov-import");
  var fileElemOvImport = el("file-elem-ov-import");

  // Buttons
  var btnExport = el("btn-export");
  var btnExportXYZ = el("btn-export-xyz");

  // GIF buttons
  var btnGifStart = el("btn-gif-start");
  var btnGifStop = el("btn-gif-stop");
  var btnGifDownload = el("btn-gif-download");
  var lblGifStatus = el("lbl-gif-status");

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
    center_A: [0, 0, 0], // anchor in Å (XY; Z kept for completeness)
    canvas_w: 400,
    canvas_h: 400,
  };

  const PAN_STEP_PX = 20;

  function getRotateStepDeg() {
    return 5;
  }

  var DEFAULT_UI = {
    mode: "TEM",
    zoom: 0.1,
    contrast: 1.35,
    background_gray: 127,
    blur_sigma: 1.0,
    noise_enabled: true,
    noise_stddev: 3.0,
    invert: false,
    draw_bonds_flag: false,
    hide_front: false,
    focus_z: 0.5,
    dof_strength: 0.0,
    scanlines: false,
    tip_sharpness: 0.5,
    clip_lo: "",
    clip_hi: "",
    bond_width: 6,
    bond_amplitude: 0.4,
    canvas_w: 400,
    canvas_h: 400,
    tip_highlight: false,
    tip_project: false,
    scale_bar: true,
  };

  // UI microscopy mode preset (does NOT change renderer math in Wave 1)
  var mode = "TEM";

  // Provider-based state (lazy-ready)
  var provider = null;
  var providerMeta = null;
  var title = "";

  // Scene state: objects + active selection.
  var scene = {
    objects: [],
    activeId: null,
  };

  function makeSceneObjectId() {
    return (
      "scene_" +
      Date.now().toString(36) +
      "_" +
      Math.random().toString(36).slice(2, 8)
    );
  }

  function normalizeClipZValue(v, fallback) {
    var n = Number.isFinite(v) ? v : parseFloat(v);
    if (!Number.isFinite(n)) n = fallback;
    if (!Number.isFinite(n)) n = 0;
    if (n < 0) n = 0;
    if (n > 1) n = 1;
    return n;
  }

  function normalizeClipZ(clip) {
    var src = clip && typeof clip === "object" ? clip : {};
    var out = {
      enabled: !!src.enabled,
      min: normalizeClipZValue(src.min, 0),
      max: normalizeClipZValue(src.max, 1),
      mode: "norm",
    };
    if (out.min > out.max) {
      var t = out.min;
      out.min = out.max;
      out.max = t;
    }
    return out;
  }

  function defaultObjectEdits() {
    return {
      posA: { x: 0, y: 0, z: 0 },
      rotDeg: { x: 0, y: 0, z: 0 },
      clipZ: { enabled: false, min: 0.0, max: 1.0, mode: "norm" },
    };
  }

  function cloneObjectEdits(edits) {
    var src = edits && typeof edits === "object" ? edits : {};
    var pos = src.posA && typeof src.posA === "object" ? src.posA : {};
    var rot = src.rotDeg && typeof src.rotDeg === "object" ? src.rotDeg : {};
    return {
      posA: {
        x: Number.isFinite(pos.x) ? pos.x : parseFloat(pos.x) || 0,
        y: Number.isFinite(pos.y) ? pos.y : parseFloat(pos.y) || 0,
        z: Number.isFinite(pos.z) ? pos.z : parseFloat(pos.z) || 0,
      },
      rotDeg: {
        x: Number.isFinite(rot.x) ? rot.x : parseFloat(rot.x) || 0,
        y: Number.isFinite(rot.y) ? rot.y : parseFloat(rot.y) || 0,
        z: Number.isFinite(rot.z) ? rot.z : parseFloat(rot.z) || 0,
      },
      clipZ: normalizeClipZ(src.clipZ),
    };
  }

  function ensureSceneObjectEdits(obj) {
    if (!obj) return defaultObjectEdits();
    obj.edits = cloneObjectEdits(obj.edits);
    return obj.edits;
  }

  function getActiveSceneObject() {
    var id = scene && scene.activeId;
    if (!id || !scene || !Array.isArray(scene.objects)) return null;
    for (var i = 0; i < scene.objects.length; i++) {
      if (scene.objects[i] && scene.objects[i].id === id)
        return scene.objects[i];
    }
    return null;
  }

  function findSceneObjectIndexById(id) {
    if (!scene || !Array.isArray(scene.objects)) return -1;
    for (var i = 0; i < scene.objects.length; i++) {
      if (scene.objects[i] && scene.objects[i].id === id) return i;
    }
    return -1;
  }

  function updateActiveSourceIndicator() {
    // Scene-driven in C1/C2.
  }

  // SMILES error de-dup (avoid alert spam for same input)
  var lastSmilesErrorToken = "";

  // SMILES atom identity cache (fixes cases where SMILES provider returns wrong/empty Z/sym).
  // Filled on successful SMILES build; applied before rendering + element overrides table.
  var _smilesAtomZCache = { token: "", Zs: null, syms: null };

  // Active file is materialized as a temporary blob URL only for the active scene object.
  var file_state = {
    kind: null,
    name: "",
    url: null,
    isBlob: false,
  };

  // ---- Element overrides (global per-element multipliers) ----
  // { "H": { size: 1.0, dark: 1.0 }, "O": { size: 1.2, dark: 0.8 } }
  var elementOverrides = {};
  var _elemList = []; // currently displayed elements (symbols)
  var _elemListKey = ""; // cache key for quick equality
  var _elemReprAtoms = Object.create(null); // sym -> representative atom (from current view)
  var _elemPreviewCanvasBySym = Object.create(null); // sym -> preview canvas
  var _elemPreviewRafPending = false;
  var _elemPreviewQueueAll = false;
  var _elemPreviewQueueSyms = Object.create(null);

  // ---- GIF recording ----
  // Recording is session-based: Start begins (or resumes) a session, Stop pauses it,
  // Download finalizes and unlocks size controls.
  var gifRec = new GifRecorder(cvs, ctx, {
    maxColors: 256,
    repeat: 0,
    maxFrames: 900,
  });

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
        lblGifStatus.textContent = "";
      } else {
        var n = gifRec.getFrameCount();
        lblGifStatus.textContent = recording
          ? tr("gifRecPrefix", "GIF: recording… frames=") + n
          : tr("gifPausePrefix", "GIF: paused. frames=") +
            n +
            tr("gifPauseSuffix", " (Download to finish)");
      }
    }
  }

  // i18n helper (uses globals set by translations loader in index.html/help.html)
  function tr(key, fallback) {
    try {
      var d = window.__I18N_DICT;
      if (d && typeof d[key] === "string") return d[key];
    } catch (e) {}
    return fallback;
  }

  function getSmilesInvalidMsg() {
    return tr("smiles_invalid", "Некоректний SMILES / Invalid SMILES");
  }

  function setSmilesErrorVisible(on, msg) {
    if (!lblSmilesErr) return;
    try {
      if (on && typeof msg === "string") lblSmilesErr.textContent = msg;
      // handle both mechanisms: class-based hiding (ui_hidden) and inline style
      try {
        lblSmilesErr.classList.toggle("ui_hidden", !on);
      } catch (e1) {}
      try {
        lblSmilesErr.style.display = on ? "" : "none";
      } catch (e2) {}
    } catch (e) {}
  }

  function maybeAlertSmilesError(msg, token, explicit) {
    if (!explicit) return;
    var t = String(token || "").trim();
    if (t && t === lastSmilesErrorToken) return;
    try {
      window.alert(msg);
    } catch (e) {}
    if (t) lastSmilesErrorToken = t;
  }

  // Validate typed SMILES using RDKit (prevents silent fallback to 'O' inside phases/builders).
  // Returns true if SMILES looks valid (or RDKit not available), false if invalid.
  async function rdkitValidateSmiles(smiles) {
    smiles = trimStr(smiles);
    if (!smiles) return true; // empty is allowed -> will fall back to 'O' by design
    var mol = null;
    try {
      // RDKit is initialized in index.html as window.RDKitReady = initRDKitModule(...)
      if (!window.RDKitReady) return true;
      var rdkit = await window.RDKitReady;
      if (!rdkit) return true;

      // rdkit.js API: get_mol(smiles) throws on invalid
      if (typeof rdkit.get_mol === "function") mol = rdkit.get_mol(smiles);
      else if (typeof rdkit.getMol === "function") mol = rdkit.getMol(smiles);
      else return true; // unknown API -> don't block

      if (!mol) return false;

      var nAtoms = 1;
      try {
        if (typeof mol.get_num_atoms === "function")
          nAtoms = mol.get_num_atoms();
        else if (typeof mol.getNumAtoms === "function")
          nAtoms = mol.getNumAtoms();
      } catch (e1) {
        nAtoms = 1;
      }

      try {
        mol.delete();
      } catch (e2) {}
      return (nAtoms | 0) > 0;
    } catch (e) {
      try {
        if (mol) mol.delete();
      } catch (e3) {}
      return false;
    }
  }

  function _normElemSym(s) {
    if (typeof s !== "string") return null;
    s = s.trim();
    if (!s) return null;
    if (s.length === 1) return s.toUpperCase();
    return s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
  }

  function _rdkitParseAtomSymsFromMolblock(molblock) {
    if (typeof molblock !== "string" || !molblock) return null;
    var lines = molblock.replace(/\r/g, "").split("\n");
    if (!lines.length) return null;

    var hasV3000 = false;
    for (var i = 0; i < Math.min(lines.length, 12); i++) {
      if (String(lines[i] || "").indexOf("V3000") >= 0) {
        hasV3000 = true;
        break;
      }
    }

    if (hasV3000) {
      var out3 = [];
      var inAtoms = false;
      for (var j = 0; j < lines.length; j++) {
        var ln3 = String(lines[j] || "").trim();
        if (!ln3) continue;
        if (ln3.indexOf("M  V30 BEGIN ATOM") === 0) {
          inAtoms = true;
          continue;
        }
        if (ln3.indexOf("M  V30 END ATOM") === 0) break;
        if (!inAtoms) continue;
        if (ln3.indexOf("M  V30 ") !== 0) continue;
        var parts3 = ln3.split(/\s+/);
        if (parts3.length < 5) continue;
        var sym3 = _normElemSym(parts3[3]);
        if (sym3) out3.push(sym3);
      }
      return out3.length ? out3 : null;
    }

    if (lines.length < 4) return null;
    var counts = String(lines[3] || "");
    var atomCount = parseInt(counts.slice(0, 3), 10);
    if (!Number.isFinite(atomCount) || atomCount <= 0) {
      var m = counts.match(/^\s*(\d+)/);
      atomCount = m ? parseInt(m[1], 10) : 0;
    }
    atomCount = atomCount | 0;
    if (atomCount <= 0) return null;
    if (lines.length < 4 + atomCount) return null;

    var out2 = new Array(atomCount);
    for (var k = 0; k < atomCount; k++) {
      var ln2 = String(lines[4 + k] || "");
      var sym2 = _normElemSym(ln2.slice(31, 34));
      if (!sym2) {
        var parts2 = ln2.trim().split(/\s+/);
        if (parts2.length >= 4) sym2 = _normElemSym(parts2[3]);
      }
      out2[k] = sym2 || null;
    }
    return out2;
  }

  function _rdkitParseAtomInfoFromJson(jsonText) {
    if (typeof jsonText !== "string" || !jsonText) return null;

    var root = null;
    try {
      root = JSON.parse(jsonText);
    } catch (e0) {
      return null;
    }
    if (!root || typeof root !== "object") return null;

    var mol = root;
    if (Array.isArray(root.molecules) && root.molecules.length)
      mol = root.molecules[0];
    if (!mol || typeof mol !== "object") return null;

    var atoms = Array.isArray(mol.atoms)
      ? mol.atoms
      : mol.mol && Array.isArray(mol.mol.atoms)
        ? mol.mol.atoms
        : null;
    if (!Array.isArray(atoms) || !atoms.length) return null;

    var syms = new Array(atoms.length);
    var zs = new Array(atoms.length);
    var any = false;

    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i] || {};
      var z = 0;
      if (Number.isFinite(a.z)) z = a.z | 0;
      else if (Number.isFinite(a.atomicNum)) z = a.atomicNum | 0;
      else if (Number.isFinite(a.atomic_num)) z = a.atomic_num | 0;

      var sym = null;
      if (z > 0 && _Z2SYM && z < _Z2SYM.length) sym = _Z2SYM[z];
      if (!sym) {
        sym = _normElemSym(
          a.symbol || a.sym || a.el || a.element || a.atom || a.label,
        );
      }
      if (!z && sym && _SYM2Z && _SYM2Z[sym]) z = _SYM2Z[sym] | 0;

      syms[i] = sym || null;
      zs[i] = z > 0 ? z : 0;
      if (syms[i] || zs[i]) any = true;
    }

    return any ? { Zs: zs, syms: syms } : null;
  }

  // Extract per-atom atomic numbers/symbols from RDKit for SMILES.
  // Returns { Zs:Array<int>|null, syms:Array<string>|null } or null.
  async function rdkitGetAtomZs(smiles) {
    smiles = trimStr(smiles);
    if (!smiles) smiles = "O";

    var mol = null;
    try {
      if (!window.RDKitReady) return null;
      var rdkit = await window.RDKitReady;
      if (!rdkit) return null;

      if (typeof rdkit.get_mol === "function") mol = rdkit.get_mol(smiles);
      else if (typeof rdkit.getMol === "function") mol = rdkit.getMol(smiles);
      else return null;

      if (!mol) return null;

      var n = 0;
      try {
        if (typeof mol.get_num_atoms === "function") n = mol.get_num_atoms();
        else if (typeof mol.getNumAtoms === "function") n = mol.getNumAtoms();
      } catch (e0) {
        n = 0;
      }

      n = n | 0;
      if (n <= 0) {
        try {
          mol.delete();
        } catch (e1) {}
        return null;
      }

      var zs = new Array(n);
      var syms = new Array(n);
      var okCount = 0;

      for (var i = 0; i < n; i++) {
        var at = null;
        try {
          if (typeof mol.get_atom_with_idx === "function")
            at = mol.get_atom_with_idx(i);
          else if (typeof mol.getAtomWithIdx === "function")
            at = mol.getAtomWithIdx(i);
        } catch (e2) {
          at = null;
        }

        var z = 0;
        var sym = null;
        if (at) {
          try {
            if (typeof at.get_atomic_num === "function")
              z = at.get_atomic_num();
            else if (typeof at.getAtomicNum === "function")
              z = at.getAtomicNum();
          } catch (e3) {
            z = 0;
          }
          try {
            if (typeof at.get_symbol === "function") sym = at.get_symbol();
            else if (typeof at.getSymbol === "function") sym = at.getSymbol();
          } catch (e4) {
            sym = null;
          }
          try {
            if (typeof at.delete === "function") at.delete();
          } catch (e5) {}
        }

        z = (z | 0) > 0 ? z | 0 : 0;
        sym = _normElemSym(sym);
        if (!sym && z > 0 && _Z2SYM && z < _Z2SYM.length) sym = _Z2SYM[z];
        if (!z && sym && _SYM2Z && _SYM2Z[sym]) z = _SYM2Z[sym] | 0;

        zs[i] = z;
        syms[i] = sym || null;
        if (z > 0 || sym) okCount++;
      }

      if (okCount < n) {
        var infoJson = null;
        try {
          if (typeof mol.get_json === "function") {
            infoJson = _rdkitParseAtomInfoFromJson(mol.get_json());
          } else if (typeof mol.getJson === "function") {
            infoJson = _rdkitParseAtomInfoFromJson(mol.getJson());
          }
        } catch (e6) {
          infoJson = null;
        }
        if (
          infoJson &&
          Array.isArray(infoJson.syms) &&
          infoJson.syms.length === n
        ) {
          for (var j = 0; j < n; j++) {
            if (!syms[j] && infoJson.syms[j]) syms[j] = infoJson.syms[j];
            if (
              (!zs[j] || zs[j] <= 0) &&
              infoJson.Zs &&
              (infoJson.Zs[j] | 0) > 0
            ) {
              zs[j] = infoJson.Zs[j] | 0;
            }
          }
        }
      }

      var needMolblock = false;
      for (var t = 0; t < n; t++) {
        if (!syms[t]) {
          needMolblock = true;
          break;
        }
      }
      if (needMolblock) {
        var molSyms = null;
        try {
          if (typeof mol.get_molblock === "function")
            molSyms = _rdkitParseAtomSymsFromMolblock(mol.get_molblock());
          else if (typeof mol.getMolblock === "function")
            molSyms = _rdkitParseAtomSymsFromMolblock(mol.getMolblock());
        } catch (e7) {
          molSyms = null;
        }
        if (Array.isArray(molSyms) && molSyms.length === n) {
          for (var q = 0; q < n; q++) {
            if (!syms[q] && molSyms[q]) syms[q] = molSyms[q];
          }
        }
      }

      for (var u = 0; u < n; u++) {
        if ((!zs[u] || zs[u] <= 0) && syms[u] && _SYM2Z && _SYM2Z[syms[u]]) {
          zs[u] = _SYM2Z[syms[u]] | 0;
        }
        if (!syms[u] && (zs[u] | 0) > 0 && _Z2SYM && zs[u] < _Z2SYM.length) {
          syms[u] = _Z2SYM[zs[u] | 0];
        }
      }

      try {
        mol.delete();
      } catch (e8) {}
      return { Zs: zs, syms: syms };
    } catch (e) {
      try {
        if (mol) mol.delete();
      } catch (e9) {}
      return null;
    }
  }

  function _smilesSourceIdx(a) {
    if (!a || typeof a !== "object") return -1;
    var keys = [
      "atom_idx",
      "atomIdx",
      "atomIndex",
      "rdkit_idx",
      "rdkitIdx",
      "source_idx",
      "sourceIdx",
      "idx",
      "index",
    ];
    for (var i = 0; i < keys.length; i++) {
      var v = a[keys[i]];
      if (!Number.isFinite(v)) v = parseInt(v, 10);
      if (Number.isFinite(v)) return v | 0;
    }
    return -1;
  }

  function _applySmilesAtomIdentity(a, z, sym) {
    if (!a) return;

    var zi = Number.isFinite(z) ? z | 0 : parseInt(z, 10);
    if (!Number.isFinite(zi)) zi = 0;
    zi = zi | 0;

    var s = _normElemSym(sym);
    if (!s && zi > 0 && _Z2SYM && zi < _Z2SYM.length) s = _Z2SYM[zi];
    if ((!zi || zi <= 0) && s && _SYM2Z && _SYM2Z[s]) zi = _SYM2Z[s] | 0;

    if (zi > 0) a.Z = zi;
    if (s) a.sym = s;
  }

  function applySmilesAtomZsIfAvailable(atoms) {
    if (!atoms || !atoms.length) return;
    if (!_smilesAtomZCache) return;

    var hasZs =
      Array.isArray(_smilesAtomZCache.Zs) && _smilesAtomZCache.Zs.length > 0;
    var hasSyms =
      Array.isArray(_smilesAtomZCache.syms) &&
      _smilesAtomZCache.syms.length > 0;
    if (!hasZs && !hasSyms) return;

    var tok = trimStr(tbSmiles ? tbSmiles.value : "") || "O";
    if (_smilesAtomZCache.token && tok !== _smilesAtomZCache.token) return;

    var zs = hasZs ? _smilesAtomZCache.Zs : null;
    var syms = hasSyms ? _smilesAtomZCache.syms : null;
    var srcLen = Math.max(zs ? zs.length : 0, syms ? syms.length : 0) | 0;
    if (srcLen <= 0) return;

    var mapByIdx = true;
    var used = Object.create(null);
    var mapped = 0;
    for (var i = 0; i < atoms.length; i++) {
      var srcIdx = _smilesSourceIdx(atoms[i]);
      if (srcIdx < 0 || srcIdx >= srcLen || used[srcIdx]) {
        mapByIdx = false;
        break;
      }
      used[srcIdx] = true;
      _applySmilesAtomIdentity(
        atoms[i],
        zs ? zs[srcIdx] : 0,
        syms ? syms[srcIdx] : null,
      );
      mapped++;
    }
    if (mapByIdx && mapped === atoms.length) return;

    if (atoms.length === srcLen) {
      for (var j = 0; j < atoms.length; j++) {
        _applySmilesAtomIdentity(
          atoms[j],
          zs ? zs[j] : 0,
          syms ? syms[j] : null,
        );
      }
      return;
    }

    var srcNonH = [];
    for (var k = 0; k < srcLen; k++) {
      var sym = syms && syms[k] ? _normElemSym(syms[k]) : null;
      var z = zs ? zs[k] | 0 : 0;
      if (!sym && z > 0 && _Z2SYM && z < _Z2SYM.length) sym = _Z2SYM[z];
      if (!z && sym && _SYM2Z && _SYM2Z[sym]) z = _SYM2Z[sym] | 0;
      if (sym === "H" || z === 1) continue;
      srcNonH.push(k);
    }

    if (srcNonH.length === atoms.length) {
      for (var q = 0; q < atoms.length; q++) {
        var src = srcNonH[q];
        _applySmilesAtomIdentity(
          atoms[q],
          zs ? zs[src] : 0,
          syms ? syms[src] : null,
        );
      }
    }
  }

  function splitSmilesCandidates(text) {
    var src = String(text || "").replace(/\r/g, "");
    var out = [];
    var cur = "";
    var depthRound = 0;
    var depthSquare = 0;
    var depthCurly = 0;

    for (var i = 0; i < src.length; i++) {
      var ch = src.charAt(i);

      if (ch === "(") depthRound++;
      else if (ch === ")") depthRound = Math.max(0, depthRound - 1);
      else if (ch === "[") depthSquare++;
      else if (ch === "]") depthSquare = Math.max(0, depthSquare - 1);
      else if (ch === "{") depthCurly++;
      else if (ch === "}") depthCurly = Math.max(0, depthCurly - 1);

      var topLevel = depthRound === 0 && depthSquare === 0 && depthCurly === 0;
      var isSeparator = false;

      if (topLevel) {
        if (
          ch === "\n" ||
          ch === "," ||
          ch === ";" ||
          ch === "|" ||
          ch === "/" ||
          ch === "\\"
        ) {
          isSeparator = true;
        } else if (/\s/.test(ch)) {
          isSeparator = true;
        }
      }

      if (isSeparator) {
        var token = trimStr(cur);
        if (token) out.push(token);
        cur = "";
      } else {
        cur += ch;
      }
    }

    var last = trimStr(cur);
    if (last) out.push(last);
    return out;
  }

  async function addSmilesObjectFromInput(opts) {
    opts = opts || {};
    var explicit = opts.explicit === true;
    var raw = trimStr(tbSmiles ? tbSmiles.value : "");
    if (!raw) return;

    var vals = splitSmilesCandidates(raw);
    if (!vals.length) return;

    var invalid = [];
    var added = [];
    for (var i = 0; i < vals.length; i++) {
      var val = vals[i];
      var valid = await rdkitValidateSmiles(val);
      if (!valid) {
        invalid.push(val);
        continue;
      }
      var obj = createSceneObjectFromSmiles(val);
      if (!obj) {
        invalid.push(val);
        continue;
      }
      scene.objects.push(obj);
      added.push(obj);
    }

    if (invalid.length) {
      var msg =
        invalid.length === 1
          ? getSmilesInvalidMsg() + ": " + invalid[0]
          : getSmilesInvalidMsg() + " (" + invalid.length + ")";
      setSmilesErrorVisible(true, msg);
      maybeAlertSmilesError(msg, invalid.join(" | "), explicit);
    } else {
      setSmilesErrorVisible(false);
      lastSmilesErrorToken = "";
    }

    if (!added.length) return;

    renderSceneList();
    saveSceneMetaToSession();
    activateSceneObject(added[added.length - 1].id, {
      render: true,
      rebuildCfg: { explicit: explicit },
    });
  }

  // ---- utils ----
  function trimStr(s) {
    return (s || "").trim();
  }

  // Element symbol helper.
  // Many providers keep only Z, so we use a minimal Z->symbol table as a fallback.
  var _Z2SYM = [
    null,
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
  ];

  // Fast symbol -> Z lookup (best effort). Used for preview when providers omit Z.
  var _SYM2Z = (function () {
    var m = Object.create(null);
    for (var z = 1; z < _Z2SYM.length; z++) {
      var s = _Z2SYM[z];
      if (s && typeof s === "string") m[s] = z;
    }
    return m;
  })();

  function atomSymbol(a) {
    if (!a) return null;

    // Prefer atomic number -> symbol (SMILES/RDKit sometimes mis-populates .sym)
    var Z = a.Z;
    if (!Number.isFinite(Z)) Z = parseInt(Z, 10);
    if (Number.isFinite(Z)) {
      Z = Z | 0;
      if (Z > 0 && Z < _Z2SYM.length) {
        var zs = _Z2SYM[Z];
        if (zs) return zs;
      }
    }

    var s =
      a.sym ||
      a.el ||
      a.symbol ||
      a.element ||
      a.Sym ||
      a.El ||
      a.Symbol ||
      a.Element;
    if (typeof s === "string") {
      s = s.trim();
      if (s) return s;
    }
    return null;
  }

  function exportAtomSymbol(a) {
    var sym = atomSymbol(a);
    sym = _normElemSym(sym);
    if (sym) return sym;
    var z = a && Number.isFinite(a.Z) ? a.Z | 0 : parseInt(a && a.Z, 10);
    if (Number.isFinite(z) && z > 0 && z < _Z2SYM.length && _Z2SYM[z])
      return _Z2SYM[z];
    return "X";
  }

  function clampMul(x, lo, hi) {
    var v = typeof x === "number" ? x : parseFloat(x);
    if (!Number.isFinite(v)) return 1.0;
    if (lo != null) v = Math.max(lo, v);
    if (hi != null) v = Math.min(hi, v);
    return v;
  }

  function gatherElementsFromAtoms(atoms) {
    var mapZ = Object.create(null);
    var mapAtom = Object.create(null);
    var out = [];
    if (!atoms || !atoms.length) return out;
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i];
      var sym = atomSymbol(a);
      if (!sym) continue;
      if (mapZ[sym] != null) continue;
      var z = a && Number.isFinite(a.Z) ? a.Z | 0 : null;
      mapZ[sym] = z;
      mapAtom[sym] = a;
      out.push(sym);
    }

    _elemReprAtoms = mapAtom;
    out.sort(function (a, b) {
      var za = mapZ[a];
      var zb = mapZ[b];
      if (Number.isFinite(za) && Number.isFinite(zb) && za !== zb)
        return za - zb;
      if (Number.isFinite(za) && !Number.isFinite(zb)) return -1;
      if (!Number.isFinite(za) && Number.isFinite(zb)) return 1;
      return String(a).localeCompare(String(b));
    });
    return out;
  }

  // ---- Element preview (single-atom) ----
  function scheduleElemPreview(sym) {
    if (!sym) return;
    if (!_elemPreviewCanvasBySym || !_elemPreviewCanvasBySym[sym]) return;
    _elemPreviewQueueSyms[sym] = true;
    _scheduleElemPreviewRAF();
  }

  function scheduleElemPreviewAll() {
    _elemPreviewQueueAll = true;
    _scheduleElemPreviewRAF();
  }

  function _scheduleElemPreviewRAF() {
    if (!elemOvDetails || !elemOvDetails.open) return;
    if (_elemPreviewRafPending) return;
    _elemPreviewRafPending = true;
    requestAnimationFrame(function () {
      _elemPreviewRafPending = false;
      if (_elemPreviewQueueAll) {
        _elemPreviewQueueAll = false;
        _elemPreviewQueueSyms = Object.create(null);
        renderAllElementPreviews();
        return;
      }

      var syms = [];
      for (var k in _elemPreviewQueueSyms) {
        if (Object.prototype.hasOwnProperty.call(_elemPreviewQueueSyms, k))
          syms.push(k);
      }
      _elemPreviewQueueSyms = Object.create(null);
      for (var i = 0; i < syms.length; i++) renderElementPreview(syms[i]);
    });
  }

  function renderAllElementPreviews() {
    if (!elemOvDetails || !elemOvDetails.open) return;
    if (!_elemList || !_elemList.length) return;
    for (var i = 0; i < _elemList.length; i++)
      renderElementPreview(_elemList[i]);
  }

  function renderElementPreview(sym) {
    var cv = _elemPreviewCanvasBySym ? _elemPreviewCanvasBySym[sym] : null;
    if (!cv) return;
    var pw = cv.width | 0;
    var ph = cv.height | 0;
    if (pw <= 0 || ph <= 0) return;

    var pctx = null;
    try {
      pctx = cv.getContext("2d", { willReadFrequently: true });
    } catch (e) {
      pctx = null;
    }
    if (!pctx) return;

    // Representative Z (best effort)
    var repr = _elemReprAtoms ? _elemReprAtoms[sym] : null;
    var Z = repr && Number.isFinite(repr.Z) ? repr.Z | 0 : _SYM2Z[sym] || 6;
    if (!Number.isFinite(Z) || Z <= 0) Z = 6;

    // Preview zoom: keep within a visible band
    var aPerPx = rngZoom ? parseFloat(rngZoom.value) : view.angstroms_per_pixel;
    if (!Number.isFinite(aPerPx)) aPerPx = 0.1;
    aPerPx = Math.max(0.03, Math.min(0.25, aPerPx));

    var cam = make_view_state({
      img_size: [ph, pw],
      angstroms_per_pixel: aPerPx,
      pan_px: [0, 0],
      rotZ_rad: 0,
      center_A: [0, 0, 0],
      center_mode: "bbox",
    });
    if (cam) {
      cam.canvas_w = pw;
      cam.canvas_h = ph;
      cam.rot_deg = 0;
      cam.rotX_rad = 0;
      cam.rotY_rad = 0;
      cam.tilt_x_deg = 0;
      cam.tilt_y_deg = 0;
    }

    var atom = { x: 0, y: 0, z: 0, Z: Z, sym: sym };

    var opts = {
      bonds: null,
      img_size: [ph, pw],
      angstroms_per_pixel: aPerPx,
      blur_sigma: rngBlur ? parseFloat(rngBlur.value) : 1.0,
      background_gray: getBackgroundGray(),
      invert: cbInvert ? !!cbInvert.checked : false,
      noise_stddev:
        cbNoise && cbNoise.checked && rngNoise
          ? parseFloat(rngNoise.value)
          : 0.0,
      contrast: rngContrast ? parseFloat(rngContrast.value) : 1.0,
      compose_mode: "sum",
      draw_bonds_flag: false,
      camera: cam,
      bond_wave_width_px: rngBwidth ? parseFloat(rngBwidth.value) : 6,
      bond_wave_amplitude: rngBamp ? parseFloat(rngBamp.value) : 0.4,
      low_clip:
        tbClipLo && tbClipLo.value !== "" ? parseFloat(tbClipLo.value) : null,
      high_clip:
        tbClipHi && tbClipHi.value !== "" ? parseFloat(tbClipHi.value) : null,
      focal_z: 0,
      dof_strength: mode === "TEM" && rngDof ? parseFloat(rngDof.value) : 0.0,
      hide_front: false,
      show_scale_bar: false,
      canvasCtx: pctx,

      // UI-only
      mode: mode,
      scanlines: cbScanlines ? !!cbScanlines.checked : false,
      tip_sharpness: rngTip ? parseFloat(rngTip.value) : 0.5,

      // Wave 3
      element_overrides: elementOverrides,
    };

    try {
      render_image([atom], opts);
    } catch (e2) {
      // silent
    }
  }

  function rebuildElementOverridesTable(elements) {
    if (!elemOvTableWrap) return;

    elements = elements || [];
    _elemList = elements.slice(0);

    // Clear
    while (elemOvTableWrap.firstChild)
      elemOvTableWrap.removeChild(elemOvTableWrap.firstChild);

    _elemPreviewCanvasBySym = Object.create(null);

    // Build table
    var table = document.createElement("table");
    var thead = document.createElement("thead");
    var hr = document.createElement("tr");

    function th(txt) {
      var h = document.createElement("th");
      h.textContent = txt;
      return h;
    }

    hr.appendChild(th(tr("ui.elemOverrides.col.element", "Element")));
    hr.appendChild(th(tr("ui.elemOverrides.col.size", "Size ×")));
    hr.appendChild(th(tr("ui.elemOverrides.col.dark", "Darkness ×")));
    hr.appendChild(th(tr("ui.elemOverrides.col.preview", "Preview")));
    thead.appendChild(hr);
    table.appendChild(thead);

    var tbody = document.createElement("tbody");

    function makeNumInput(sym, field) {
      var inp = document.createElement("input");
      inp.type = "number";
      inp.step = "0.1";
      inp.min = "0.2";

      var cur =
        elementOverrides && elementOverrides[sym]
          ? elementOverrides[sym]
          : null;
      var v = 1.0;
      if (cur && field === "size" && Number.isFinite(cur.size)) v = cur.size;
      if (cur && field === "dark" && Number.isFinite(cur.dark)) v = cur.dark;
      inp.value = String(v);

      function apply() {
        // Keep minimum clamp only (no upper limit)
        var vv = clampMul(inp.value, 0.2, null);
        inp.value = String(vv);

        if (!elementOverrides) elementOverrides = {};
        if (!elementOverrides[sym])
          elementOverrides[sym] = { size: 1.0, dark: 1.0 };
        elementOverrides[sym][field] = vv;

        scheduleRender(); // rerender only (no rebuild/parsing)
        scheduleElemPreview(sym);
      }

      inp.addEventListener("input", apply);
      inp.addEventListener("change", apply);

      return inp;
    }

    for (var i = 0; i < elements.length; i++) {
      var sym = elements[i];
      var trr = document.createElement("tr");

      var td0 = document.createElement("td");
      td0.textContent = sym;
      trr.appendChild(td0);

      var td1 = document.createElement("td");
      td1.appendChild(makeNumInput(sym, "size"));
      trr.appendChild(td1);

      var td2 = document.createElement("td");
      td2.appendChild(makeNumInput(sym, "dark"));
      trr.appendChild(td2);

      var td3 = document.createElement("td");
      td3.className = "preview_cell";
      var cv = document.createElement("canvas");
      cv.width = 96;
      cv.height = 96;
      cv.className = "elem_preview_canvas";
      td3.appendChild(cv);
      _elemPreviewCanvasBySym[sym] = cv;
      trr.appendChild(td3);

      tbody.appendChild(trr);
    }

    table.appendChild(tbody);
    elemOvTableWrap.appendChild(table);

    // Previews: rerender when panel is open
    scheduleElemPreviewAll();
  }

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
    m = String(m || "")
      .trim()
      .toUpperCase();
    if (m === "TEM" || m === "STM" || m === "AFM") return m;
    return "TEM";
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
    try {
      node.classList.toggle("ui_hidden", !!hidden);
    } catch (e) {
      node.style.display = hidden ? "none" : "";
    }
  }

  function setWrapDisabled(input, disabled) {
    if (!input) return;
    input.disabled = !!disabled;

    var wrap = null;
    try {
      wrap = input.closest("label") || input.closest(".range_row");
    } catch (e) {
      wrap = null;
    }

    if (wrap) {
      try {
        wrap.classList.toggle("ui_disabled", !!disabled);
      } catch (e) {}
    } else {
      // fallback: direct style
      try {
        input.style.opacity = disabled ? "0.5" : "";
      } catch (e) {}
    }
  }

  function updateModeDescText() {
    if (!lblModeDesc) return;
    if (mode === "STM")
      lblModeDesc.textContent = tr(
        "modeDescSTM",
        "STM — поверхневий режим, focus/DoF вимкнено.",
      );
    else if (mode === "AFM")
      lblModeDesc.textContent = tr(
        "modeDescAFM",
        "AFM — поверхневий режим, focus/DoF вимкнено.",
      );
    else
      lblModeDesc.textContent = tr(
        "modeDescTEM",
        "TEM — проєкційний режим; STM/AFM — поверхневі режими, focus/DoF вимкнено.",
      );
  }

  function updateBondsLabelText() {
    if (!lblBonds) return;
    if (mode === "AFM") {
      lblBonds.textContent = tr("lblBondsAFM", "Bond contrast (AFM)");
    } else {
      lblBonds.textContent = tr("lblBonds", "Bonds");
    }
  }

  function applyModePreset(newMode, cfg) {
    cfg = cfg || {};
    var save = cfg.save !== false;
    var rerender = cfg.rerender !== false;

    mode = normMode(newMode);
    if (selMode) selMode.value = mode;

    if (save) {
      try {
        localStorage.setItem("scigentem_mode", mode);
      } catch (e) {}
    }

    // Visibility for mode-specific controls
    setHidden(rowScanlines, mode !== "STM");
    setHidden(rowTip, mode !== "AFM");

    // Global policies for this wave:
    //  - TEM: focus + DoF active
    //  - STM/AFM: focus + DoF disabled in UI; DoF is also hard-zeroed in render opts
    if (rngDof) rngDof.value = clampRangeValue(rngDof, 0.0);
    setWrapDisabled(cbHide, false);

    // Bonds label is UI-only; keep it i18n-aware
    updateBondsLabelText();

    if (mode === "TEM") {
      setWrapDisabled(rngFocus, false);
      setWrapDisabled(rngDof, false);

      if (cbBonds) cbBonds.checked = false;
      setWrapDisabled(cbBonds, true);

      if (cbInvert) cbInvert.checked = false;
      if (rngBlur) rngBlur.value = clampRangeValue(rngBlur, 1.0);

      if (cbNoise) cbNoise.checked = true;
      if (rngNoise) rngNoise.value = clampRangeValue(rngNoise, 3.0);

      if (rngContrast) rngContrast.value = clampRangeValue(rngContrast, 1.35);

      if (cbScanlines) cbScanlines.checked = false;
      setWrapDisabled(cbScanlines, true);
      setWrapDisabled(rngTip, true);
    }

    if (mode === "STM") {
      setWrapDisabled(rngFocus, true);
      setWrapDisabled(rngDof, true);

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

    if (mode === "AFM") {
      setWrapDisabled(rngFocus, true);
      setWrapDisabled(rngDof, true);

      if (cbBonds) cbBonds.checked = true;
      setWrapDisabled(cbBonds, false);

      if (cbInvert) cbInvert.checked = false;
      if (rngBlur) rngBlur.value = clampRangeValue(rngBlur, 0.8);

      if (cbNoise) cbNoise.checked = true;
      if (rngNoise) rngNoise.value = clampRangeValue(rngNoise, 1.2);

      if (rngContrast) rngContrast.value = clampRangeValue(rngContrast, 1.7);

      if (rngTip) rngTip.value = clampRangeValue(rngTip, 0.5);
      setWrapDisabled(rngTip, false);
      if (cbScanlines) cbScanlines.checked = false;
      setWrapDisabled(cbScanlines, true);
    }

    // Keep background label in sync
    if (rngBg) getBackgroundGray();

    updateModeDescText();
    if (rerender) scheduleRender();
  }

  var SETTINGS_KEY = "scigentem_settings_v1";
  var settingsSaveTimer = null;
  var lastSampleId = null;

  var SAMPLES = {
    "smiles:caffeine": {
      kind: "smiles",
      smiles: "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
    },
    "smiles:benzene": {
      kind: "smiles",
      smiles: "c1ccccc1",
    },
    "file:Graphite.cif": {
      kind: "file",
      path: new URL("../../samples/Graphite.cif", import.meta.url).href,
      name: "Graphite.cif",
    },
    "file:Corundum.poscar": {
      kind: "file",
      path: new URL("../../samples/Corundum.poscar", import.meta.url).href,
      name: "Corundum.poscar",
    },
    "file:pdb4hhb.ent": {
      kind: "file",
      path: new URL("../../samples/pdb4hhb.ent", import.meta.url).href,
      name: "pdb4hhb.ent",
    },
  };

  function applyDefaultUiToDom(cfg) {
    cfg = cfg || DEFAULT_UI;
    if (selMode) selMode.value = normMode(cfg.mode);
    if (rngZoom) rngZoom.value = String(cfg.zoom);
    if (rngContrast) rngContrast.value = String(cfg.contrast);
    if (rngBg) rngBg.value = String(cfg.background_gray);
    if (rngBlur) rngBlur.value = String(cfg.blur_sigma);
    if (cbNoise) cbNoise.checked = !!cfg.noise_enabled;
    if (rngNoise) rngNoise.value = String(cfg.noise_stddev);
    if (cbInvert) cbInvert.checked = !!cfg.invert;
    if (cbBonds) cbBonds.checked = !!cfg.draw_bonds_flag;
    if (cbHide) cbHide.checked = !!cfg.hide_front;
    if (rngFocus) rngFocus.value = String(cfg.focus_z);
    if (rngDof) rngDof.value = String(cfg.dof_strength);
    if (cbScanlines) cbScanlines.checked = !!cfg.scanlines;
    if (rngTip) rngTip.value = String(cfg.tip_sharpness);
    if (tbClipLo)
      tbClipLo.value = cfg.clip_lo == null ? "" : String(cfg.clip_lo);
    if (tbClipHi)
      tbClipHi.value = cfg.clip_hi == null ? "" : String(cfg.clip_hi);
    if (rngBwidth) rngBwidth.value = String(cfg.bond_width);
    if (rngBamp) rngBamp.value = String(cfg.bond_amplitude);
    if (cbTipHighlight) cbTipHighlight.checked = !!cfg.tip_highlight;
    if (cbTipProject) cbTipProject.checked = !!cfg.tip_project;
    if (cbScale) cbScale.checked = !!cfg.scale_bar;
    if (!(sizeLock && sizeLock.active)) {
      if (tbW) tbW.value = String(cfg.canvas_w);
      if (tbH) tbH.value = String(cfg.canvas_h);
    }
  }

  function resetUiSettings() {
    clearSettingsStorage();
    applyDefaultUiToDom(DEFAULT_UI);
    applyModePreset(DEFAULT_UI.mode, { save: false, rerender: false });
    getBackgroundGray();
    saveSettingsDebounced();
    scheduleRender();
  }

  function resetAllObjectEdits() {
    if (!scene || !Array.isArray(scene.objects) || !scene.objects.length)
      return;
    for (var i = 0; i < scene.objects.length; i++) {
      if (!scene.objects[i]) continue;
      scene.objects[i].edits = defaultObjectEdits();
    }
    updateActiveEditsUi();
    renderSceneList();
    saveSceneMetaToSession();
    scheduleRender();
  }

  function collectSettingsFromUI() {
    var cfg = {};
    cfg.mode = mode;
    cfg.activeSource = scene && scene.activeId ? "scene" : "smiles";

    if (rngZoom) cfg.zoom = parseFloat(rngZoom.value);
    if (rngContrast) cfg.contrast = parseFloat(rngContrast.value);
    if (rngBg) cfg.bg = parseFloat(rngBg.value);
    if (rngBlur) cfg.blur = parseFloat(rngBlur.value);
    if (rngNoise) cfg.noiseStd = parseFloat(rngNoise.value);
    if (rngDof) cfg.dof = parseFloat(rngDof.value);
    if (rngBwidth) cfg.bondWidth = parseFloat(rngBwidth.value);
    if (rngBamp) cfg.bondAmp = parseFloat(rngBamp.value);
    // NOTE: bond length min/max was removed (not needed)

    if (cbNoise) cfg.noiseEnabled = !!cbNoise.checked;
    if (cbInvert) cfg.invert = !!cbInvert.checked;
    if (cbScale) cfg.scaleBar = !!cbScale.checked;
    if (cbScanlines) cfg.scanlines = !!cbScanlines.checked;
    if (cbTipHighlight) cfg.tipHighlight = !!cbTipHighlight.checked;
    if (cbTipProject) cfg.tipProject = !!cbTipProject.checked;
    if (cbBonds) cfg.bonds = !!cbBonds.checked;
    if (cbHide) cfg.hideFront = !!cbHide.checked;

    if (rngFocus) cfg.focus = parseFloat(rngFocus.value);

    cfg.tiles = {
      nx: tbNx ? parseInt(tbNx.value, 10) || 1 : 1,
      ny: tbNy ? parseInt(tbNy.value, 10) || 1 : 1,
      nz: tbNz ? parseInt(tbNz.value, 10) || 1 : 1,
    };

    if (tbW) cfg.canvasW = parseInt(tbW.value, 10) || 0;
    if (tbH) cfg.canvasH = parseInt(tbH.value, 10) || 0;

    if (tbSmiles) cfg.smilesText = tbSmiles.value || "";

    cfg.lastSampleId = lastSampleId;

    return cfg;
  }

  function applySettingsFromObject(cfg) {
    if (!cfg || typeof cfg !== "object") return;

    var targetMode = cfg.mode ? normMode(cfg.mode) : null;
    if (targetMode) {
      if (selMode) selMode.value = targetMode;
      applyModePreset(targetMode, { save: false, rerender: false });
    }

    function setNum(input, val) {
      if (!input || typeof val === "undefined" || val === null) return;
      var num = typeof val === "number" ? val : parseFloat(val);
      if (!isFinite(num)) return;
      input.value = String(num);
    }

    function setBool(cb, val) {
      if (!cb || typeof val !== "boolean") return;
      cb.checked = val;
    }

    setNum(rngZoom, cfg.zoom);
    setNum(rngContrast, cfg.contrast);
    setNum(rngBg, cfg.bg);
    setNum(rngBlur, cfg.blur);
    setNum(rngNoise, cfg.noiseStd);
    setNum(rngDof, cfg.dof);
    setNum(rngBwidth, cfg.bondWidth);
    setNum(rngBamp, cfg.bondAmp);

    setBool(cbNoise, cfg.noiseEnabled);
    setBool(cbInvert, cfg.invert);
    setBool(cbScale, cfg.scaleBar);
    setBool(cbScanlines, cfg.scanlines);
    setBool(cbTipHighlight, cfg.tipHighlight);
    setBool(cbTipProject, cfg.tipProject);
    setBool(cbBonds, cfg.bonds);
    setBool(cbHide, cfg.hideFront);

    setNum(rngFocus, cfg.focus);

    if (cfg.tiles) {
      setNum(tbNx, cfg.tiles.nx);
      setNum(tbNy, cfg.tiles.ny);
      setNum(tbNz, cfg.tiles.nz);
    }

    if (typeof cfg.canvasW !== "undefined") setNum(tbW, cfg.canvasW);
    if (typeof cfg.canvasH !== "undefined") setNum(tbH, cfg.canvasH);

    if (typeof cfg.smilesText === "string" && tbSmiles) {
      tbSmiles.value = cfg.smilesText;
    }

    if (cfg.activeSource === "smiles") {
      updateActiveSourceIndicator();
    }

    if (cfg.lastSampleId) {
      lastSampleId = cfg.lastSampleId;
    }
  }

  function saveSettingsDebounced() {
    if (!window.localStorage) return;

    var payload;
    try {
      payload = JSON.stringify(collectSettingsFromUI());
    } catch (e) {
      return;
    }

    if (settingsSaveTimer !== null) {
      clearTimeout(settingsSaveTimer);
    }

    settingsSaveTimer = setTimeout(function () {
      settingsSaveTimer = null;
      try {
        window.localStorage.setItem(SETTINGS_KEY, payload);
      } catch (e) {
        // Ignore quota and access errors
      }
    }, 250);
  }

  function restoreSettingsFromStorage() {
    if (!window.localStorage) return false;

    var raw;
    try {
      raw = window.localStorage.getItem(SETTINGS_KEY);
    } catch (e) {
      return false;
    }

    if (!raw) return false;

    var cfg;
    try {
      cfg = JSON.parse(raw);
    } catch (e) {
      return false;
    }

    if (!cfg || typeof cfg !== "object") return false;

    applySettingsFromObject(cfg);
    return true;
  }

  function clearSettingsStorage() {
    if (!window.localStorage) return;
    try {
      window.localStorage.removeItem(SETTINGS_KEY);
    } catch (e) {
      // Ignore
    }
  }

  function setSamplesStatus(msgKey) {
    if (!lblSamplesStatus) return;
    if (!msgKey) {
      lblSamplesStatus.textContent = "";
      return;
    }
    var msg = tr(msgKey) || "";
    lblSamplesStatus.textContent = msg;
  }

  async function applySampleById(sampleId) {
    if (!sampleId || !SAMPLES[sampleId]) return;

    var sample = SAMPLES[sampleId];
    lastSampleId = sampleId;

    if (sample.kind === "smiles") {
      if (!tbSmiles) return;
      tbSmiles.value = sample.smiles || "";
      updateActiveSourceIndicator();
      setSmilesErrorVisible(false);
      addSmilesObjectFromInput({ explicit: true });
      setSamplesStatus(null);
      if (selSamples) selSamples.value = "__placeholder__";
      saveSettingsDebounced();
      return;
    }

    if (sample.kind === "file" && sample.path && window.fetch) {
      setSamplesStatus("sampleLoading");
      try {
        var resp = await fetch(sample.path);
        if (!resp.ok) {
          throw new Error("HTTP " + resp.status);
        }
        var text = await resp.text();
        var fileName =
          sample.name || sample.path.split("/").pop() || "sample.dat";
        var obj = createSceneObjectFromFile(fileName, text);
        if (obj) addSceneObject(obj, { makeActive: true, render: true });
      } catch (e) {
        console.error("Failed to load sample", e);
        alert(
          tr("sampleFetchFailed") ||
            "Sample files not found. Put them into /samples/.",
        );
      } finally {
        setSamplesStatus(null);
        if (selSamples) selSamples.value = "__placeholder__";
      }
      saveSettingsDebounced();
    }
  }

  function clearCanvas(bg) {
    ctx.save();
    ctx.fillStyle = "rgb(" + bg + "," + bg + "," + bg + ")";
    ctx.fillRect(0, 0, cvs.width, cvs.height);
    ctx.restore();
  }

  function setTitle(t) {
    title = t || "";
    if (lblTitle) lblTitle.textContent = "Simulated EM — " + title;
  }

  function fileKindFromName(name) {
    var low = String(name || "").toLowerCase();
    // tolerate things like "file.cif;..." or "file.cif?x" in the #fragment
    if (low.indexOf(".cif") >= 0) return "cif";
    if (low.indexOf(".csv") >= 0) return "csv";
    if (low.indexOf(".xyz") >= 0) return "xyz";
    if (low.indexOf(".json") >= 0) return "json";
    if (low.indexOf(".mol") >= 0) return "mol";
    if (low.indexOf(".poscar") >= 0) return "poscar";
    if (low.indexOf(".contcar") >= 0) return "poscar";
    if (low.indexOf(".vasp") >= 0) return "poscar";
    // common VASP names without extension
    if (low === "poscar" || low === "contcar") return "poscar";
    return null;
  }

  function revokeFileUrl() {
    if (file_state.isBlob && file_state.url) {
      try {
        URL.revokeObjectURL(file_state.url);
      } catch (e) {}
    }
    file_state.kind = null;
    file_state.name = "";
    file_state.url = null;
    file_state.isBlob = false;
  }

  function parseTilesSize() {
    function readTileN(tb) {
      var v = tb ? parseInt(tb.value, 10) : 1;
      if (!Number.isFinite(v) || v <= 0) v = 1;
      return v | 0;
    }
    var nx = readTileN(tbNx);
    var ny = readTileN(tbNy);
    var nz = readTileN(tbNz);

    // normalize UI (so empty/invalid snaps back to 1)
    if (tbNx) tbNx.value = String(nx);
    if (tbNy) tbNy.value = String(ny);
    if (tbNz) tbNz.value = String(nz);

    return String(nx) + "x" + String(ny) + "x" + String(nz);
  }

  function applyPeriodicSizeToSpec(spec, sizeStr) {
    if (!sizeStr) return spec;

    var s = String(spec || "").trim();
    if (!s) return s;

    var semi = s.indexOf(";");
    if (semi < 0) return s + ";size=" + sizeStr;

    var base = s.slice(0, semi).trim();
    var rest = s.slice(semi + 1).trim();
    var items = rest ? rest.split(";") : [];

    var kept = [];
    for (var i = 0; i < items.length; i++) {
      var it = trimStr(items[i]);
      if (!it) continue;
      var itLow = it.toLowerCase();
      if (itLow.indexOf("size=") === 0) continue;
      kept.push(it);
    }
    kept.push("size=" + sizeStr);

    return base + ";" + kept.join(";");
  }

  function getSceneMetaForSession() {
    var out = [];
    if (!scene || !Array.isArray(scene.objects)) return out;
    for (var i = 0; i < scene.objects.length; i++) {
      var o = scene.objects[i];
      if (!o) continue;
      var edits = cloneObjectEdits(o.edits);
      out.push({
        id: o.id,
        name: o.name || "",
        visible: o.visible !== false,
        createdAt: o.createdAt || 0,
        kind: o.source && o.source.kind ? o.source.kind : "",
        format: o.source && o.source.format ? o.source.format : "",
        edits: edits,
      });
    }
    return out;
  }

  function saveSceneMetaToSession() {
    try {
      sessionStorage.setItem(
        "scigentem_scene_meta_v2",
        JSON.stringify({
          activeId: scene.activeId,
          objects: getSceneMetaForSession(),
        }),
      );
    } catch (e) {}
  }

  function updateTilesUiForActiveObject() {
    var active = getActiveSceneObject();
    var tiles =
      active && Array.isArray(active.tiles) ? active.tiles : [1, 1, 1];
    if (tbNx) tbNx.value = String((tiles[0] | 0) > 0 ? tiles[0] | 0 : 1);
    if (tbNy) tbNy.value = String((tiles[1] | 0) > 0 ? tiles[1] | 0 : 1);
    if (tbNz) tbNz.value = String((tiles[2] | 0) > 0 ? tiles[2] | 0 : 1);
  }

  function fmtShortNum(v) {
    var n = Number(v);
    if (!Number.isFinite(n)) n = 0;
    var s = n.toFixed(1);
    if (s === "-0.0") s = "0.0";
    return s;
  }

  function updateActiveEditsUi() {
    var active = getActiveSceneObject();
    var edits = active ? ensureSceneObjectEdits(active) : defaultObjectEdits();
    var clipZ = edits.clipZ || normalizeClipZ();
    if (lblMoveValues) {
      lblMoveValues.textContent =
        "x=" +
        fmtShortNum(edits.posA.x) +
        ", y=" +
        fmtShortNum(edits.posA.y) +
        ", z=" +
        fmtShortNum(edits.posA.z) +
        " Å";
    }
    if (lblRotateValues) {
      lblRotateValues.textContent =
        "x=" +
        fmtShortNum(edits.rotDeg.x) +
        ", y=" +
        fmtShortNum(edits.rotDeg.y) +
        ", z=" +
        fmtShortNum(edits.rotDeg.z) +
        "°";
    }
    if (cbClipZEnabled) cbClipZEnabled.checked = !!clipZ.enabled;
    if (rngClipZMin) rngClipZMin.value = clipZ.min.toFixed(2);
    if (rngClipZMax) rngClipZMax.value = clipZ.max.toFixed(2);
    if (lblClipZValues) {
      lblClipZValues.textContent = tr("ui.clipZ.values", "min={min} max={max}")
        .replace("{min}", clipZ.min.toFixed(2))
        .replace("{max}", clipZ.max.toFixed(2));
    }
    var disabled = !active;
    [
      btnMoveUp,
      btnMoveDown,
      btnMoveLeft,
      btnMoveRight,
      btnMoveHome,
      btnMoveZp,
      btnMoveZm,
      btnRotateXp,
      btnRotateXm,
      btnRotateYp,
      btnRotateYm,
      btnRotateHome,
      btnRotateZp,
      btnRotateZm,
      btnResetObjectEdits,
      cbClipZEnabled,
      rngClipZMin,
      rngClipZMax,
      btnClipZReset,
    ].forEach(function (x) {
      if (x) x.disabled = disabled;
    });
    if (rngClipZMin) rngClipZMin.disabled = disabled || !clipZ.enabled;
    if (rngClipZMax) rngClipZMax.disabled = disabled || !clipZ.enabled;
  }

  function getSceneObjectLabel(obj) {
    if (!obj) return "";
    var kind = obj.source && obj.source.kind ? obj.source.kind : "";
    if (kind === "smiles") return "SMILES";
    var fmt =
      obj.source && obj.source.format
        ? String(obj.source.format).toUpperCase()
        : "";
    return fmt;
  }

  function renderSceneList() {
    if (!sceneList || !sceneEmpty) return;
    sceneList.innerHTML = "";
    var list = scene && Array.isArray(scene.objects) ? scene.objects : [];
    sceneEmpty.textContent = tr("ui.scene.empty", "No objects");
    sceneEmpty.style.display = list.length ? "none" : "";
    if (!list.length) {
      updateActiveEditsUi();
      return;
    }
    for (var i = 0; i < list.length; i++) {
      (function (obj) {
        var row = document.createElement("div");
        row.className =
          "scene_item" + (obj.id === scene.activeId ? " is_active" : "");

        var vis = document.createElement("input");
        vis.type = "checkbox";
        vis.checked = obj.visible !== false;
        vis.title = tr("ui.scene.visible", "Visible");
        vis.addEventListener("click", function (ev) {
          ev.stopPropagation();
        });
        vis.addEventListener("change", function (ev) {
          obj.visible = !!vis.checked;
          saveSceneMetaToSession();
          updateSceneStatsUi();
          scheduleRender();
          ev.stopPropagation();
        });
        row.appendChild(vis);

        var main = document.createElement("div");
        main.className = "scene_item_main";
        main.addEventListener("click", function () {
          activateSceneObject(obj.id);
        });
        var name = document.createElement("div");
        name.className = "scene_item_name";
        name.textContent = obj.name || "Untitled";
        var meta = document.createElement("div");
        meta.className = "scene_item_meta";
        meta.textContent = getSceneObjectLabel(obj);
        main.appendChild(name);
        main.appendChild(meta);
        row.appendChild(main);

        var btn = document.createElement("button");
        btn.type = "button";
        btn.className = "scene_remove_btn";
        btn.textContent = "✕";
        btn.title = tr("ui.scene.remove", "Remove");
        btn.addEventListener("click", function (ev) {
          ev.stopPropagation();
          removeSceneObject(obj.id);
        });
        row.appendChild(btn);
        sceneList.appendChild(row);
      })(list[i]);
    }
    updateSceneStatsUi();
    updateActiveEditsUi();
  }

  function updateSceneStatsUi() {
    if (!sceneStats) return;
    var list = scene && Array.isArray(scene.objects) ? scene.objects : [];
    var total = list.length;
    var rendered = 0;
    for (var i = 0; i < list.length; i++) {
      if (list[i] && list[i].visible !== false) rendered++;
    }
    sceneStats.textContent = tr(
      "ui.scene.stats",
      "Rendered: {rendered} / {total}",
    )
      .replace("{rendered}", String(rendered))
      .replace("{total}", String(total));
  }

  function activateSceneObject(id, cfg) {
    cfg = cfg || {};
    if (!id) {
      scene.activeId = null;
      renderSceneList();
      updateTilesUiForActiveObject();
      saveSceneMetaToSession();
      if (cfg.render !== false) scheduleRender();
      return;
    }
    scene.activeId = id;
    updateTilesUiForActiveObject();
    renderSceneList();
    saveSceneMetaToSession();
    if (cfg.render !== false) rebuildAndRender(cfg.rebuildCfg || {});
  }

  function addSceneObject(obj, cfg) {
    cfg = cfg || {};
    ensureSceneObjectEdits(obj);
    scene.objects.push(obj);
    renderSceneList();
    saveSceneMetaToSession();
    if (cfg.makeActive !== false) {
      activateSceneObject(obj.id, {
        render: cfg.render !== false,
        rebuildCfg: cfg.rebuildCfg || {},
      });
    }
  }

  function removeSceneObject(id) {
    var idx = findSceneObjectIndexById(id);
    if (idx < 0) return;
    var wasActive = scene.activeId === id;
    scene.objects.splice(idx, 1);
    if (!wasActive) {
      renderSceneList();
      saveSceneMetaToSession();
      return;
    }
    if (!scene.objects.length) {
      activateSceneObject(null, { render: true });
      return;
    }
    var nextIdx = Math.min(idx - 1, scene.objects.length - 1);
    if (nextIdx < 0) nextIdx = scene.objects.length - 1;
    activateSceneObject(scene.objects[nextIdx].id, { render: true });
  }

  function createSceneObjectFromFile(name, text) {
    var kind = fileKindFromName(name);
    if (!kind) return null;
    return {
      id: makeSceneObjectId(),
      name: name || "file",
      source: { kind: "file", format: kind, text: text || "" },
      tiles: isPeriodicFileKind(kind, name) ? [1, 1, 1] : null,
      visible: true,
      createdAt: Date.now(),
      edits: defaultObjectEdits(),
    };
  }

  function createSceneObjectFromSmiles(smiles) {
    var val = trimStr(smiles);
    if (!val) return null;
    return {
      id: makeSceneObjectId(),
      name: "SMILES: " + val,
      source: { kind: "smiles", smiles: val },
      tiles: null,
      visible: true,
      createdAt: Date.now(),
      edits: defaultObjectEdits(),
    };
  }

  async function handleLocalFiles(files) {
    if (!files || !files.length) return;
    var added = [];
    for (var i = 0; i < files.length; i++) {
      var file = files[i];
      if (!file) continue;
      var kind = fileKindFromName(file.name);
      if (!kind) continue;
      var text = "";
      try {
        text = await file.text();
      } catch (e) {
        continue;
      }
      var obj = createSceneObjectFromFile(file.name, text);
      if (obj) {
        scene.objects.push(obj);
        added.push(obj);
      }
    }
    if (!added.length) return;
    renderSceneList();
    saveSceneMetaToSession();
    activateSceneObject(added[added.length - 1].id, { render: true });
  }

  function resetActiveObjectPosition() {
    var active = getActiveSceneObject();
    if (!active) return;
    var edits = ensureSceneObjectEdits(active);
    edits.posA.x = 0;
    edits.posA.y = 0;
    edits.posA.z = 0;
    updateActiveEditsUi();
    saveSceneMetaToSession();
  }

  function resetActiveObjectRotation() {
    var active = getActiveSceneObject();
    if (!active) return;
    var edits = ensureSceneObjectEdits(active);
    edits.rotDeg.x = 0;
    edits.rotDeg.y = 0;
    edits.rotDeg.z = 0;
    updateActiveEditsUi();
    saveSceneMetaToSession();
  }

  function resetActiveObjectEdits() {
    var active = getActiveSceneObject();
    if (!active) return;
    active.edits = defaultObjectEdits();
    updateActiveEditsUi();
    saveSceneMetaToSession();
  }

  function resetActiveObjectClipZ() {
    var active = getActiveSceneObject();
    if (!active) return;
    var edits = ensureSceneObjectEdits(active);
    edits.clipZ = normalizeClipZ({ enabled: false, min: 0, max: 1 });
    updateActiveEditsUi();
    saveSceneMetaToSession();
  }

  function filterAtomsAndBondsByClipZ(atoms, bonds, clipZ) {
    var clip = normalizeClipZ(clipZ);
    if (!clip.enabled || !atoms || !atoms.length) {
      return { atoms: atoms || [], bonds: bonds };
    }
    var zmin = Infinity;
    var zmax = -Infinity;
    for (var i = 0; i < atoms.length; i++) {
      var zi = atoms[i] && atoms[i].z != null ? Number(atoms[i].z) || 0 : 0;
      if (zi < zmin) zmin = zi;
      if (zi > zmax) zmax = zi;
    }
    if (!Number.isFinite(zmin) || !Number.isFinite(zmax)) {
      return { atoms: atoms || [], bonds: Array.isArray(bonds) ? [] : bonds };
    }
    var span = zmax - zmin;
    var lo = span <= 1e-12 ? zmin : zmin + clip.min * span;
    var hi = span <= 1e-12 ? zmax : zmin + clip.max * span;
    var outAtoms = [];
    var indexMap = new Array(atoms.length);
    for (var j = 0; j < atoms.length; j++) {
      var a = atoms[j];
      var z = a && a.z != null ? Number(a.z) || 0 : 0;
      if (z < lo || z > hi) {
        indexMap[j] = -1;
        continue;
      }
      indexMap[j] = outAtoms.length;
      outAtoms.push(a);
    }
    var outBonds = bonds;
    if (Array.isArray(bonds)) {
      outBonds = [];
      for (var k = 0; k < bonds.length; k++) {
        var b = bonds[k];
        if (!b || b.length < 2) continue;
        var i0 = b[0] | 0;
        var i1 = b[1] | 0;
        var ni0 = indexMap[i0];
        var ni1 = indexMap[i1];
        if (ni0 == null || ni0 < 0 || ni1 == null || ni1 < 0) continue;
        outBonds.push([ni0, ni1, b[2] || 1]);
      }
      try {
        if (bonds && bonds.__tem_guessed) outBonds.__tem_guessed = true;
      } catch (e) {}
    }
    return { atoms: outAtoms, bonds: outBonds };
  }

  function mutateActiveObjectEdits(mutator) {
    var active = getActiveSceneObject();
    if (!active) return;
    var edits = ensureSceneObjectEdits(active);
    if (typeof mutator === "function") mutator(edits);
    updateActiveEditsUi();
    saveSceneMetaToSession();
    scheduleRender();
  }

  function bboxCenterFromAtoms(atoms) {
    if (!atoms || !atoms.length) return [0, 0, 0];
    var xmin = Infinity,
      xmax = -Infinity,
      ymin = Infinity,
      ymax = -Infinity,
      zmin = Infinity,
      zmax = -Infinity;
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i] || {};
      var x = Number(a.x) || 0;
      var y = Number(a.y) || 0;
      var z = a.z != null ? Number(a.z) || 0 : 0;
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
    }
    return [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)];
  }

  function applySceneObjectEditsToAtoms(obj, atoms) {
    if (!obj || !atoms || !atoms.length) return atoms;
    var edits = ensureSceneObjectEdits(obj);
    var tx = Number(edits.posA.x) || 0;
    var ty = Number(edits.posA.y) || 0;
    var tz = Number(edits.posA.z) || 0;
    var rx = ((Number(edits.rotDeg.x) || 0) * Math.PI) / 180.0;
    var ry = ((Number(edits.rotDeg.y) || 0) * Math.PI) / 180.0;
    var rz = ((Number(edits.rotDeg.z) || 0) * Math.PI) / 180.0;
    if (
      Math.abs(tx) < 1e-12 &&
      Math.abs(ty) < 1e-12 &&
      Math.abs(tz) < 1e-12 &&
      Math.abs(rx) < 1e-12 &&
      Math.abs(ry) < 1e-12 &&
      Math.abs(rz) < 1e-12
    ) {
      return atoms;
    }
    var center = bboxCenterFromAtoms(atoms);
    var cx = center[0],
      cy = center[1],
      cz = center[2];
    var cosx = Math.cos(rx),
      sinx = Math.sin(rx);
    var cosy = Math.cos(ry),
      siny = Math.sin(ry);
    var cosz = Math.cos(rz),
      sinz = Math.sin(rz);
    var out = new Array(atoms.length);
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i] || {};
      var x = (Number(a.x) || 0) - cx;
      var y = (Number(a.y) || 0) - cy;
      var z = (a.z != null ? Number(a.z) || 0 : 0) - cz;

      if (Math.abs(rx) > 1e-12) {
        var y1 = y * cosx - z * sinx;
        var z1 = y * sinx + z * cosx;
        y = y1;
        z = z1;
      }
      if (Math.abs(ry) > 1e-12) {
        var x2 = x * cosy + z * siny;
        var z2 = -x * siny + z * cosy;
        x = x2;
        z = z2;
      }
      if (Math.abs(rz) > 1e-12) {
        var x3 = x * cosz - y * sinz;
        var y3 = x * sinz + y * cosz;
        x = x3;
        y = y3;
      }
      out[i] = Object.assign({}, a, {
        x: cx + x + tx,
        y: cy + y + ty,
        z: cz + z + tz,
      });
    }
    return out;
  }

  function applyActiveObjectEditsToAtoms(atoms) {
    return applySceneObjectEditsToAtoms(getActiveSceneObject(), atoms);
  }

  function makeSceneObjectBuildKey(obj) {
    if (!obj || !obj.source) return "";
    var base = obj.id || "";
    if (obj.source.kind === "smiles") {
      return base + "|smiles|" + String(obj.source.smiles || "");
    }
    var tiles = Array.isArray(obj.tiles) ? obj.tiles.join("x") : "";
    return (
      base +
      "|file|" +
      String(obj.name || "") +
      "|fmt=" +
      String(obj.source.format || "") +
      "|tiles=" +
      tiles +
      "|len=" +
      String((obj.source.text || "").length)
    );
  }

  async function buildSceneObjectProvider(obj, cfg) {
    cfg = cfg || {};
    if (!obj || !obj.source) return false;
    var key = makeSceneObjectBuildKey(obj);
    if (obj._provider && obj._providerMeta && obj._providerKey === key)
      return true;

    var spec = "";
    var kind = null;
    var isSmiles = false;
    var smilesToken = "";
    var blobUrl = null;
    try {
      if (obj.source.kind === "file") {
        kind = obj.source.format;
        blobUrl = URL.createObjectURL(
          new Blob([String(obj.source.text || "")], { type: "text/plain" }),
        );
        spec =
          blobUrl + "#" + (obj.name || (kind ? "scene." + kind : "scene.dat"));
        if (isPeriodicFileKind(kind, obj.name)) {
          var tilesArr = Array.isArray(obj.tiles) ? obj.tiles : [1, 1, 1];
          var tiles = [
            Math.max(1, parseInt(tilesArr[0], 10) || 1),
            Math.max(1, parseInt(tilesArr[1], 10) || 1),
            Math.max(1, parseInt(tilesArr[2], 10) || 1),
          ].join("x");
          spec = applyPeriodicSizeToSpec(spec, tiles);
        }
      } else {
        smilesToken = obj.source ? trimStr(obj.source.smiles) : "";
        spec = smilesToken || "O";
        kind = "smiles";
        isSmiles = true;
      }

      if (isSmiles) {
        var _smOk = true;
        if (smilesToken) {
          _smOk = await rdkitValidateSmiles(smilesToken);
          if (!_smOk) return false;
        }
      }

      var newProvider = await build_system_from_input(spec);
      var newMeta =
        newProvider && typeof newProvider.getMeta === "function"
          ? newProvider.getMeta()
          : null;
      obj._provider = newProvider;
      obj._providerMeta = newMeta;
      obj._providerKey = key;

      if (isSmiles) {
        try {
          var _info = await rdkitGetAtomZs(spec);
          obj._smilesAtomZCache = {
            token: spec,
            Zs: _info && Array.isArray(_info.Zs) ? _info.Zs : null,
            syms: _info && Array.isArray(_info.syms) ? _info.syms : null,
          };
        } catch (e1) {
          obj._smilesAtomZCache = { token: spec, Zs: null, syms: null };
        }
      }
      return true;
    } finally {
      if (blobUrl) {
        try {
          URL.revokeObjectURL(blobUrl);
        } catch (e2) {}
      }
    }
  }

  function applyStoredSmilesIdentity(obj, atoms) {
    if (!obj || !atoms || !atoms.length) return;
    var prev = _smilesAtomZCache;
    try {
      if (obj._smilesAtomZCache) _smilesAtomZCache = obj._smilesAtomZCache;
      applySmilesAtomZsIfAvailable(atoms);
    } catch (e) {}
    _smilesAtomZCache = prev;
  }

  function makeBaseRenderViewState(w, h) {
    view.canvas_w = w;
    view.canvas_h = h;
    view.img_size = [h, w];
    view.angstroms_per_pixel = rngZoom
      ? parseFloat(rngZoom.value)
      : view.angstroms_per_pixel;
    var viewState = make_view_state({
      img_size: view.img_size,
      angstroms_per_pixel: view.angstroms_per_pixel,
      pan_px: view.pan_px,
      rotZ_rad: (view.rot_deg * Math.PI) / 180.0,
      center_A: view.center_A || [0, 0, 0],
      center_mode: "bbox",
    });
    viewState.canvas_w = view.canvas_w;
    viewState.canvas_h = view.canvas_h;
    viewState.rot_deg = view.rot_deg;
    viewState.rotX_rad = (view.tilt_x_deg * Math.PI) / 180.0;
    viewState.rotY_rad = (view.tilt_y_deg * Math.PI) / 180.0;
    viewState.tilt_x_deg = view.tilt_x_deg;
    viewState.tilt_y_deg = view.tilt_y_deg;
    return viewState;
  }

  function cloneUsedCamera(cam, w, h, fallbackCenterA) {
    var usedCamera = Object.assign({}, cam || {});
    if (
      (!usedCamera.center_A || !Array.isArray(usedCamera.center_A)) &&
      Array.isArray(fallbackCenterA)
    )
      usedCamera.center_A = fallbackCenterA;
    usedCamera.canvas_w = w;
    usedCamera.canvas_h = h;
    usedCamera.rot_deg = view.rot_deg;
    usedCamera.rotX_rad = (view.tilt_x_deg * Math.PI) / 180.0;
    usedCamera.rotY_rad = (view.tilt_y_deg * Math.PI) / 180.0;
    usedCamera.tilt_x_deg = view.tilt_x_deg;
    usedCamera.tilt_y_deg = view.tilt_y_deg;
    return usedCamera;
  }

  function computeFocalZForAtoms(atoms, usedCamera) {
    if (!atoms || !atoms.length) return 0;
    var focus_norm = rngFocus ? parseFloat(rngFocus.value) : 0.5;
    if (!Number.isFinite(focus_norm)) focus_norm = 0.5;
    var cxA =
      usedCamera && usedCamera.center_A ? usedCamera.center_A[0] || 0 : 0;
    var cyA =
      usedCamera && usedCamera.center_A ? usedCamera.center_A[1] || 0 : 0;
    var czA =
      usedCamera && usedCamera.center_A ? usedCamera.center_A[2] || 0 : 0;
    var ax =
      usedCamera && Number.isFinite(usedCamera.rotX_rad)
        ? usedCamera.rotX_rad
        : 0;
    var ay =
      usedCamera && Number.isFinite(usedCamera.rotY_rad)
        ? usedCamera.rotY_rad
        : 0;
    function zAfterModelTilt(a) {
      var x = (a.x || 0) - cxA;
      var y = (a.y || 0) - cyA;
      var z = (a.z != null ? a.z : 0) - czA;
      if (ax) {
        var c = Math.cos(ax),
          s = Math.sin(ax);
        var y1 = y * c - z * s;
        var z1 = y * s + z * c;
        y = y1;
        z = z1;
      }
      if (ay) {
        var c2 = Math.cos(ay),
          s2 = Math.sin(ay);
        var x2 = x * c2 + z * s2;
        var z2 = -x * s2 + z * c2;
        x = x2;
        z = z2;
      }
      return z + czA;
    }
    var useTilt = Math.abs(ax) > 1e-12 || Math.abs(ay) > 1e-12;
    var zmin = useTilt
      ? zAfterModelTilt(atoms[0])
      : atoms[0].z != null
        ? atoms[0].z
        : 0;
    var zmax = zmin;
    for (var i = 1; i < atoms.length; i++) {
      var zc = useTilt
        ? zAfterModelTilt(atoms[i])
        : atoms[i].z != null
          ? atoms[i].z
          : 0;
      if (zc < zmin) zmin = zc;
      if (zc > zmax) zmax = zc;
    }
    return zmin + focus_norm * (zmax - zmin);
  }

  function makePerObjectRenderOpts(w, h, atoms, bonds, usedCamera) {
    return {
      bonds: bonds,
      img_size: [h, w],
      angstroms_per_pixel: view.angstroms_per_pixel,
      blur_sigma: rngBlur ? parseFloat(rngBlur.value) : 1.0,
      background_gray: getBackgroundGray(),
      invert: false,
      noise_stddev: 0.0,
      contrast: 1.0,
      compose_mode: "sum",
      draw_bonds_flag: cbBonds ? !!cbBonds.checked : true,
      camera: usedCamera,
      bond_wave_width_px: rngBwidth ? parseFloat(rngBwidth.value) : 6,
      bond_wave_amplitude: rngBamp ? parseFloat(rngBamp.value) : 0.4,
      low_clip: null,
      high_clip: null,
      focal_z: computeFocalZForAtoms(atoms, usedCamera),
      dof_strength: mode === "TEM" && rngDof ? parseFloat(rngDof.value) : 0.0,
      hide_front: cbHide ? !!cbHide.checked : false,
      show_scale_bar: false,
      scale_bar_corner: "bl",
      scale_bar_margin_px: 12,
      canvasCtx: null,
      mode: mode,
      scanlines: cbScanlines ? !!cbScanlines.checked : false,
      tip_sharpness: rngTip ? parseFloat(rngTip.value) : 0.5,
      element_overrides: elementOverrides,
    };
  }

  function composeGrayInto(target, src) {
    if (!target || !src || target.length !== src.length) return;
    for (var i = 0; i < target.length; i++) {
      if (src[i] < target[i]) target[i] = src[i];
    }
  }

  var _randnSpare = null;
  function randnMain() {
    if (_randnSpare != null) {
      var v = _randnSpare;
      _randnSpare = null;
      return v;
    }
    var u = 0,
      v2 = 0;
    while (u === 0) u = Math.random();
    while (v2 === 0) v2 = Math.random();
    var mag = Math.sqrt(-2.0 * Math.log(u));
    var ang = 2.0 * Math.PI * v2;
    _randnSpare = mag * Math.sin(ang);
    return mag * Math.cos(ang);
  }

  function applyGlobalPostToGray(buf, bg) {
    var out = new Uint8ClampedArray(buf.length);
    var contrast = rngContrast ? parseFloat(rngContrast.value) : 1.0;
    if (!Number.isFinite(contrast)) contrast = 1.0;
    var lowClip =
      el("tb-clip-lo") && el("tb-clip-lo").value !== ""
        ? parseFloat(el("tb-clip-lo").value)
        : null;
    var highClip =
      el("tb-clip-hi") && el("tb-clip-hi").value !== ""
        ? parseFloat(el("tb-clip-hi").value)
        : null;
    var addNoise = cbNoise && cbNoise.checked && rngNoise;
    var noiseStd = addNoise ? parseFloat(rngNoise.value) : 0.0;
    if (!Number.isFinite(noiseStd)) noiseStd = 0.0;
    var invert = cbInvert ? !!cbInvert.checked : false;
    for (var i = 0; i < buf.length; i++) {
      var v = Number(buf[i]);
      if (!Number.isFinite(v)) v = bg;
      v = bg + (v - bg) * contrast;
      if (lowClip != null && Number.isFinite(lowClip) && v < lowClip)
        v = lowClip;
      if (highClip != null && Number.isFinite(highClip) && v > highClip)
        v = highClip;
      if (noiseStd > 0) v += randnMain() * noiseStd;
      if (invert) v = 255 - v;
      if (v < 0) v = 0;
      else if (v > 255) v = 255;
      out[i] = v | 0;
    }
    return out;
  }

  function drawGrayFrameToCanvas(gray, w, h) {
    var imageData = ctx.createImageData(w, h);
    for (var i = 0, p = 0; i < gray.length; i++, p += 4) {
      var v = gray[i];
      imageData.data[p] = v;
      imageData.data[p + 1] = v;
      imageData.data[p + 2] = v;
      imageData.data[p + 3] = 255;
    }
    ctx.putImageData(imageData, 0, 0);
  }

  function drawSceneScaleBar(w, h) {
    if (!(cbScale && cbScale.checked)) return;
    var ap = view.angstroms_per_pixel;
    if (!(Number.isFinite(ap) && ap > 0)) return;
    var targetPx = Math.min(Math.max(w * 0.25, 80), 200);
    var targetA = targetPx * ap;
    var k = targetA > 0 ? Math.floor(Math.log10(targetA)) : 0;
    var bestA = 0,
      bestPx = 0;
    var mults = [1, 2, 5];
    for (var shift = -3; shift <= 3; shift++) {
      var base = Math.pow(10, k + shift);
      for (var mi = 0; mi < mults.length; mi++) {
        var valA = mults[mi] * base;
        var valPx = valA / ap;
        if (valPx >= 80 && valPx <= 200) {
          bestA = valA;
          bestPx = valPx;
          break;
        }
      }
      if (bestPx > 0) break;
    }
    if (!(bestPx > 0)) {
      bestPx = targetPx;
      bestA = bestPx * ap;
    }
    var x = 12;
    var y = h - 16;
    ctx.save();
    ctx.strokeStyle = "#ffffff";
    ctx.fillStyle = "#ffffff";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(x, y);
    ctx.lineTo(x + bestPx, y);
    ctx.stroke();
    ctx.font = "12px Arial, sans-serif";
    var label =
      bestA >= 10 ? Math.round(bestA) + " Å" : bestA.toFixed(1) + " Å";
    ctx.fillText(label, x, y - 6);
    ctx.restore();
  }

  async function getFinalAtomsForSceneObject(obj, w, h) {
    if (!obj || obj.visible === false) return [];
    var built = await buildSceneObjectProvider(obj, {});
    if (!built || !obj._provider || !obj._providerMeta) return [];

    var centerA = Array.isArray(obj._providerMeta.center_A)
      ? obj._providerMeta.center_A
      : view.center_A || [0, 0, 0];
    var baseViewState = makeBaseRenderViewState(w, h);
    var objViewState = Object.assign({}, baseViewState, { center_A: centerA });
    var provView = obj._provider.getView(objViewState, {
      needBonds: cbBonds ? !!cbBonds.checked : true,
    });
    var atoms = provView && provView.atomsView ? provView.atomsView : [];
    var bonds = provView && "bondsView" in provView ? provView.bondsView : null;
    if (!atoms || !atoms.length) return [];
    if (obj.source && obj.source.kind === "smiles")
      applyStoredSmilesIdentity(obj, atoms);
    atoms = applySceneObjectEditsToAtoms(obj, atoms);
    var clipped = filterAtomsAndBondsByClipZ(
      atoms,
      bonds,
      ensureSceneObjectEdits(obj).clipZ,
    );
    atoms = clipped.atoms;
    if (!atoms || !atoms.length) return [];
    return atoms;
  }

  function buildSceneExportXYZ(atoms) {
    var lines = [];
    var list = Array.isArray(atoms) ? atoms : [];
    lines.push(String(list.length));
    lines.push(
      "SciGenTEM scene export (final coords after per-object pos/rot/clip)",
    );
    for (var i = 0; i < list.length; i++) {
      var a = list[i] || {};
      var sym = exportAtomSymbol(a);
      var x = Number(a.x) || 0;
      var y = Number(a.y) || 0;
      var z = a.z != null ? Number(a.z) || 0 : 0;
      lines.push(
        sym + "  " + x.toFixed(6) + "  " + y.toFixed(6) + "  " + z.toFixed(6),
      );
    }
    return lines.join("\n") + "\n";
  }

  async function exportSceneXYZ() {
    var list = scene && Array.isArray(scene.objects) ? scene.objects : [];
    var visible = [];
    for (var i = 0; i < list.length; i++) {
      if (list[i] && list[i].visible !== false) visible.push(list[i]);
    }
    var w = clampInt(tbW ? tbW.value : 400, 128, 4096, 400);
    var h = clampInt(tbH ? tbH.value : 400, 128, 4096, 400);
    var exportAtoms = [];
    for (var j = 0; j < visible.length; j++) {
      var atoms = await getFinalAtomsForSceneObject(visible[j], w, h);
      for (var k = 0; k < atoms.length; k++) exportAtoms.push(atoms[k]);
    }
    if (!exportAtoms.length) {
      try {
        window.alert(tr("ui.export.xyz.empty", "Nothing to export"));
      } catch (e) {}
      return;
    }
    var xyz = buildSceneExportXYZ(exportAtoms);
    var blob = new Blob([xyz], { type: "chemical/x-xyz;charset=utf-8" });
    var a = document.createElement("a");
    a.download = "scigentem_scene.xyz";
    a.href = URL.createObjectURL(blob);
    a.click();
    setTimeout(function () {
      try {
        URL.revokeObjectURL(a.href);
      } catch (e) {}
    }, 1000);
  }

  async function renderSceneFrame(w, h) {
    var list = scene && Array.isArray(scene.objects) ? scene.objects : [];
    var bg = getBackgroundGray();
    var visible = [];
    for (var i = 0; i < list.length; i++) {
      if (list[i] && list[i].visible !== false) visible.push(list[i]);
    }
    updateSceneStatsUi();
    if (!visible.length) return { gray: null, allAtoms: [] };

    var globalGray = new Float32Array(w * h);
    for (var g = 0; g < globalGray.length; g++) globalGray[g] = bg;
    var allAtoms = [];
    var baseViewState = makeBaseRenderViewState(w, h);

    for (var j = 0; j < visible.length; j++) {
      var obj = visible[j];
      var built = await buildSceneObjectProvider(obj, {});
      if (!built || !obj._provider || !obj._providerMeta) continue;
      var centerA = Array.isArray(obj._providerMeta.center_A)
        ? obj._providerMeta.center_A
        : view.center_A || [0, 0, 0];
      var objViewState = Object.assign({}, baseViewState, {
        center_A: centerA,
      });
      var provView = obj._provider.getView(objViewState, {
        needBonds: cbBonds ? !!cbBonds.checked : true,
      });
      var atoms = provView && provView.atomsView ? provView.atomsView : [];
      var bonds =
        provView && "bondsView" in provView ? provView.bondsView : null;
      if (!atoms || !atoms.length) continue;
      if (obj.source && obj.source.kind === "smiles")
        applyStoredSmilesIdentity(obj, atoms);
      atoms = applySceneObjectEditsToAtoms(obj, atoms);
      var clipped = filterAtomsAndBondsByClipZ(
        atoms,
        bonds,
        ensureSceneObjectEdits(obj).clipZ,
      );
      atoms = clipped.atoms;
      bonds = clipped.bonds;
      if (!atoms || !atoms.length) continue;
      for (var ai = 0; ai < atoms.length; ai++) allAtoms.push(atoms[ai]);
      var usedCamera = cloneUsedCamera(
        provView && provView.usedCamera ? provView.usedCamera : objViewState,
        w,
        h,
        centerA,
      );
      var opts = makePerObjectRenderOpts(w, h, atoms, bonds, usedCamera);
      var frame = render_image(atoms, opts);
      if (frame && frame.length === globalGray.length)
        composeGrayInto(globalGray, frame);
    }
    return { gray: globalGray, allAtoms: allAtoms };
  }

  // ---- build + render ----
  async function buildCurrentSystem(cfg) {
    cfg = cfg || {};
    var activeObj = getActiveSceneObject();
    if (!activeObj) {
      provider = null;
      providerMeta = null;
      setTitle("");
      setSmilesErrorVisible(false);
      revokeFileUrl();
      return true;
    }
    var ok = await buildSceneObjectProvider(activeObj, cfg);
    if (!ok) return false;
    provider = activeObj._provider || null;
    providerMeta = activeObj._providerMeta || null;
    if (providerMeta && Array.isArray(providerMeta.center_A))
      view.center_A = providerMeta.center_A;
    var kind =
      activeObj.source && activeObj.source.kind === "file"
        ? activeObj.source.format || ""
        : "smiles";
    var uiTitle =
      activeObj && activeObj.name
        ? activeObj.name
        : providerMeta && providerMeta.title
          ? providerMeta.title
          : kind || "";
    setTitle(uiTitle);
    return true;
  }

  async function renderCurrent() {
    var w = clampInt(tbW ? tbW.value : 400, 128, 4096, 400);
    var h = clampInt(tbH ? tbH.value : 400, 128, 4096, 400);
    if (tbW) tbW.value = String(w);
    if (tbH) tbH.value = String(h);

    if (cvs.width !== w) cvs.width = w;
    if (cvs.height !== h) cvs.height = h;

    var sceneFrame = await renderSceneFrame(w, h);
    if (!sceneFrame || !sceneFrame.gray) {
      clearCanvas(getBackgroundGray());
      updateGifUI();
      return;
    }

    try {
      var elems = gatherElementsFromAtoms(sceneFrame.allAtoms);
      var key = elems.join("|");
      if (key !== _elemListKey) {
        _elemListKey = key;
        rebuildElementOverridesTable(elems);
      }
    } catch (e) {}

    var finalGray = applyGlobalPostToGray(sceneFrame.gray, getBackgroundGray());
    drawGrayFrameToCanvas(finalGray, w, h);
    drawSceneScaleBar(w, h);

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

    // Keep per-element previews in sync with global effects (noise/blur/contrast/etc)
    scheduleElemPreviewAll();
    saveSettingsDebounced();
  }

  var rebuildPending = false;
  var rebuildQueued = false;
  var rebuildQueuedCfg = null;
  async function rebuildAndRender(cfg) {
    // Coalesce rebuild requests instead of dropping them.
    // This prevents the "typed SMILES -> temporary invalid -> fallback O -> never updates" issue
    // when RDKit build is still running.
    if (rebuildPending) {
      rebuildQueued = true;
      if (!rebuildQueuedCfg) rebuildQueuedCfg = {};
      if (cfg && cfg.explicit) rebuildQueuedCfg.explicit = true;
      return;
    }

    rebuildPending = true;
    var ok = true;
    try {
      ok = await buildCurrentSystem(cfg);
    } catch (e) {
      // Keep UI responsive even on parse errors.
      // (Minimal error handling; no debug spam.)
      provider = null;
      providerMeta = null;
      setTitle("ERROR");
      ok = true;
    } finally {
      rebuildPending = false;
    }

    // For invalid SMILES, keep the last good image (don't rerender / don't refresh noise).
    if (ok !== false) {
      await renderCurrent();
    }

    if (rebuildQueued) {
      var qcfg = rebuildQueuedCfg || {};
      rebuildQueuedCfg = null;
      rebuildQueued = false;
      setTimeout(function () {
        rebuildAndRender(qcfg);
      }, 0);
    }
  }

  // ---- file load ----

  // ---- events ----

  // Render-only controls
  [
    rngZoom,
    rngContrast,
    rngBlur,
    rngBg,
    rngNoise,
    rngFocus,
    rngDof,
    rngBwidth,
    rngBamp,
    tbClipLo,
    tbClipHi,
    rngTip,
  ].forEach(function (x) {
    if (!x) return;
    x.oninput = function () {
      if (x === rngBg) getBackgroundGray();
      scheduleRender();
    };
  });

  [cbInvert, cbNoise, cbBonds, cbHide, cbScale, cbScanlines].forEach(
    function (x) {
      if (x) x.onchange = scheduleRender;
    },
  );

  if (selMode) {
    selMode.onchange = function () {
      applyModePreset(selMode.value, { save: true, rerender: true });
    };
  }

  if (tbW)
    tbW.onchange = function () {
      if (sizeLock && sizeLock.active) {
        tbW.value = String(sizeLock.w);
        return;
      }
      scheduleRender();
    };
  if (tbH)
    tbH.onchange = function () {
      if (sizeLock && sizeLock.active) {
        tbH.value = String(sizeLock.h);
        return;
      }
      scheduleRender();
    };

  function onTilesChanged() {
    var active = getActiveSceneObject();
    if (!active || !active.source || active.source.kind !== "file") return;
    if (!isPeriodicFileKind(active.source.format, active.name)) return;
    active.tiles = [
      Math.max(1, parseInt(tbNx && tbNx.value, 10) || 1),
      Math.max(1, parseInt(tbNy && tbNy.value, 10) || 1),
      Math.max(1, parseInt(tbNz && tbNz.value, 10) || 1),
    ];
    saveSceneMetaToSession();
    rebuildAndRender();
  }
  [tbNx, tbNy, tbNz].forEach(function (t) {
    if (t) t.onchange = onTilesChanged;
  });

  if (tbSmiles) {
    tbSmiles.value = trimStr(smiles_text);
    if (!trimStr(tbSmiles.value)) tbSmiles.value = "O";
    tbSmiles.oninput = function () {
      setSmilesErrorVisible(false);
      lastSmilesErrorToken = "";
    };
    tbSmiles.onchange = function () {
      setSmilesErrorVisible(false);
      lastSmilesErrorToken = "";
    };
    tbSmiles.addEventListener("keydown", function (ev) {
      if (ev.key !== "Enter") return;
      var isTextarea =
        String(tbSmiles.tagName || "").toUpperCase() === "TEXTAREA";
      if (isTextarea && ev.shiftKey) return;
      ev.preventDefault();
      addSmilesObjectFromInput({ explicit: true });
    });
  }

  if (fileInput) {
    fileInput.onchange = function () {
      var files = fileInput.files;
      if (files && files.length) handleLocalFiles(files);
      fileInput.value = "";
    };
  }

  document.addEventListener("dragover", function (ev) {
    ev.preventDefault();
  });
  document.addEventListener("drop", function (ev) {
    ev.preventDefault();
    var dt = ev.dataTransfer;
    if (!dt || !dt.files || !dt.files.length) return;
    handleLocalFiles(dt.files);
  });

  if (btnExport) {
    btnExport.onclick = function () {
      var a = document.createElement("a");
      a.download = "tem.png";
      a.href = cvs.toDataURL("image/png");
      a.click();
    };
  }

  if (btnExportXYZ) {
    btnExportXYZ.onclick = function () {
      exportSceneXYZ();
    };
  }

  // ---- Element overrides controls ----
  if (elemOvDetails) {
    elemOvDetails.addEventListener("toggle", function () {
      if (elemOvDetails.open) scheduleElemPreviewAll();
    });
  }

  function pruneOverridesForExport(ov) {
    var out = {};
    if (!ov) return out;
    for (var k in ov) {
      if (!Object.prototype.hasOwnProperty.call(ov, k)) continue;
      var it = ov[k];
      if (!it || typeof it !== "object") continue;
      var sz = Number.isFinite(it.size) ? it.size : 1.0;
      var dk = Number.isFinite(it.dark) ? it.dark : 1.0;
      // Export: omit defaults (1,1)
      if (Math.abs(sz - 1.0) < 1e-12 && Math.abs(dk - 1.0) < 1e-12) continue;
      out[k] = { size: sz, dark: dk };
    }
    return out;
  }

  function buildFullOverridesForExport() {
    var out = {};
    // Export a FULL config: include all periodic table elements.
    // If an element is not present in the current structure, we export defaults (1.0, 1.0).
    if (typeof _Z2SYM !== "undefined" && _Z2SYM && _Z2SYM.length) {
      for (var Z = 1; Z < _Z2SYM.length; Z++) {
        var sym = _Z2SYM[Z];
        if (!sym) continue;
        var it =
          elementOverrides && elementOverrides[sym]
            ? elementOverrides[sym]
            : null;
        var sz = it && Number.isFinite(it.size) ? it.size : 1.0;
        var dk = it && Number.isFinite(it.dark) ? it.dark : 1.0;
        out[sym] = { size: sz, dark: dk };
      }
    }
    // Preserve any extra keys (non-periodic pseudo-elements) if present.
    if (elementOverrides) {
      for (var k in elementOverrides) {
        if (!Object.prototype.hasOwnProperty.call(elementOverrides, k))
          continue;
        if (out[k]) continue;
        var it2 = elementOverrides[k];
        var sz2 =
          it2 && typeof it2 === "object" && Number.isFinite(it2.size)
            ? it2.size
            : 1.0;
        var dk2 =
          it2 && typeof it2 === "object" && Number.isFinite(it2.dark)
            ? it2.dark
            : 1.0;
        out[k] = { size: sz2, dark: dk2 };
      }
    }
    return out;
  }

  function exportOverridesConfig() {
    try {
      var payload = {
        version: 1,
        overrides: buildFullOverridesForExport(),
      };
      var text = JSON.stringify(payload, null, 2);
      var blob = new Blob([text], { type: "application/json" });
      var a = document.createElement("a");
      a.download = tr(
        "ui.elemOverrides.exportName",
        "scigentem_overrides.json",
      );
      a.href = URL.createObjectURL(blob);
      a.click();
      setTimeout(function () {
        try {
          URL.revokeObjectURL(a.href);
        } catch (e) {}
      }, 1000);
    } catch (e) {
      // silent
    }
  }

  function importOverridesFromObject(obj) {
    if (!obj || typeof obj !== "object") return;

    var raw =
      obj.overrides && typeof obj.overrides === "object" ? obj.overrides : null;
    if (!raw) return;

    // Merge into elementOverrides (keep unknown elements; UI shows only current ones)
    for (var sym in raw) {
      if (!Object.prototype.hasOwnProperty.call(raw, sym)) continue;
      var it = raw[sym];
      if (!it || typeof it !== "object") continue;

      var sz = 1.0;
      var dk = 1.0;

      if (it.size != null) sz = clampMul(it.size, 0.2, null);
      if (it.dark != null) dk = clampMul(it.dark, 0.2, null);

      if (!elementOverrides[sym])
        elementOverrides[sym] = { size: 1.0, dark: 1.0 };
      elementOverrides[sym].size = sz;
      elementOverrides[sym].dark = dk;
    }
  }

  async function importOverridesConfigFile(file) {
    if (!file) return false;
    var text = "";
    try {
      text = await file.text();
    } catch (e) {
      return false;
    }

    var obj = null;
    try {
      obj = JSON.parse(text);
    } catch (e2) {
      return false;
    }

    try {
      importOverridesFromObject(obj);
      // refresh UI table (for current elements) + rerender
      rebuildElementOverridesTable(_elemList || []);
      scheduleRender();
      return true;
    } catch (e3) {
      return false;
    }
  }

  if (btnElemOvReset) {
    btnElemOvReset.addEventListener("click", function () {
      elementOverrides = {};
      rebuildElementOverridesTable(_elemList || []);
      scheduleRender();
    });
  }
  if (btnElemOvExport) {
    btnElemOvExport.addEventListener("click", function () {
      exportOverridesConfig();
    });
  }
  if (btnElemOvImport) {
    btnElemOvImport.addEventListener("click", function () {
      try {
        if (fileElemOvImport) fileElemOvImport.click();
      } catch (e) {}
    });
  }
  if (fileElemOvImport) {
    fileElemOvImport.addEventListener("change", async function () {
      var f =
        fileElemOvImport.files && fileElemOvImport.files.length
          ? fileElemOvImport.files[0]
          : null;
      // reset input so same file can be re-imported
      try {
        fileElemOvImport.value = "";
      } catch (e0) {}

      if (!f) return;

      var ok = await importOverridesConfigFile(f);
      if (!ok) {
        try {
          alert(tr("ui.elemOverrides.importError", "Import error"));
        } catch (e) {}
      }
    });
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

      var a = document.createElement("a");
      a.download = "tem.gif";
      a.href = URL.createObjectURL(blob);
      a.click();

      setTimeout(function () {
        try {
          URL.revokeObjectURL(a.href);
        } catch (e) {}
      }, 1000);

      // unlock size only after download finalizes
      setSizeLock(false);
      updateGifUI();
    };
  }

  // Per-object edit controls: active object only, render-only updates.
  function getMoveStepA() {
    return 0.5;
  }

  function bindObjectEditControls() {
    if (btnMoveUp)
      btnMoveUp.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.y -= getMoveStepA();
        });
      };
    if (btnMoveDown)
      btnMoveDown.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.y += getMoveStepA();
        });
      };
    if (btnMoveLeft)
      btnMoveLeft.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.x -= getMoveStepA();
        });
      };
    if (btnMoveRight)
      btnMoveRight.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.x += getMoveStepA();
        });
      };
    if (btnMoveZp)
      btnMoveZp.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.z += getMoveStepA();
        });
      };
    if (btnMoveZm)
      btnMoveZm.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.posA.z -= getMoveStepA();
        });
      };
    if (btnMoveHome)
      btnMoveHome.onclick = function () {
        resetActiveObjectPosition();
        scheduleRender();
      };

    if (btnRotateXp)
      btnRotateXp.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.x += getRotateStepDeg();
        });
      };
    if (btnRotateXm)
      btnRotateXm.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.x -= getRotateStepDeg();
        });
      };
    if (btnRotateYm)
      btnRotateYm.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.y -= getRotateStepDeg();
        });
      };
    if (btnRotateYp)
      btnRotateYp.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.y += getRotateStepDeg();
        });
      };
    if (btnRotateZm)
      btnRotateZm.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.z -= getRotateStepDeg();
        });
      };
    if (btnRotateZp)
      btnRotateZp.onclick = function () {
        mutateActiveObjectEdits(function (edits) {
          edits.rotDeg.z += getRotateStepDeg();
        });
      };
    if (btnRotateHome)
      btnRotateHome.onclick = function () {
        resetActiveObjectRotation();
        scheduleRender();
      };
    if (btnResetObjectEdits)
      btnResetObjectEdits.onclick = function () {
        resetActiveObjectEdits();
        scheduleRender();
      };
  }

  bindObjectEditControls();

  // refresh dynamic labels on language switch
  try {
    document.addEventListener("i18n-updated", function () {
      updateGifUI();
      updateModeDescText();
      updateBondsLabelText();
      renderSceneList();
      updateActiveEditsUi();
      try {
        rebuildElementOverridesTable(_elemList || []);
      } catch (e) {}
    });
  } catch (e) {}

  // ---- init ----
  var restored = false;
  try {
    restored = restoreSettingsFromStorage();
  } catch (e) {
    restored = false;
  }

  if (!restored) {
    try {
      var storedMode = normMode(localStorage.getItem("scigentem_mode"));
      applyModePreset(storedMode, { save: false, rerender: false });
    } catch (e) {
      applyModePreset("TEM", { save: false, rerender: false });
    }
  }

  // Ensure Samples dropdown is reset to placeholder
  if (selSamples && selSamples.value !== "__placeholder__") {
    selSamples.value = "__placeholder__";
  }

  // Attach handlers that depend on helpers defined above
  if (selSamples) {
    selSamples.addEventListener("change", function () {
      var v = selSamples.value;
      if (!v || v === "__placeholder__") return;
      applySampleById(v);
    });
  }

  if (btnResetSettings) {
    btnResetSettings.addEventListener("click", function () {
      resetUiSettings();
    });
  }

  if (btnResetAllObjectEdits) {
    btnResetAllObjectEdits.addEventListener("click", function () {
      resetAllObjectEdits();
    });
  }

  getBackgroundGray();
  updateGifUI();
  renderSceneList();
  updateActiveEditsUi();

  if (tbSmiles && trimStr(tbSmiles.value)) {
    await addSmilesObjectFromInput({ explicit: false });
  } else {
    scheduleRender();
  }
}
