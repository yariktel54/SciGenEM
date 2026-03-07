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
import { make_view_state } from "./camera.js";
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

  // Element overrides (per-element multipliers)
  var elemOvTableWrap = el("elem-override-table");
  var elemOvDetails = el("elem-override-details");
  var btnElemOvReset = el("btn-elem-ov-reset");
  var btnElemOvExport = el("btn-elem-ov-export");
  var btnElemOvImport = el("btn-elem-ov-import");
  var fileElemOvImport = el("file-elem-ov-import");

  // Buttons
  var btnExport = el("btn-export");

  // GIF buttons
  var btnGifStart = el("btn-gif-start");
  var btnGifStop = el("btn-gif-stop");
  var btnGifDownload = el("btn-gif-download");
  var lblGifStatus = el("lbl-gif-status");

  // View (stage) buttons
  var btnViewUp = el("btn-view-up");
  var btnViewDown = el("btn-view-down");
  var btnViewLeft = el("btn-view-left");
  var btnViewRight = el("btn-view-right");
  var btnViewRotL = el("btn-view-rotl");
  var btnViewRotR = el("btn-view-rotr");
  var btnViewReset = el("btn-view-reset");
  var btnModelX90 = el("btn-model-x90");
  var btnModelY90 = el("btn-model-y90");

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
  const ROT_STEP_DEG = 5;
  const TILT_STEP_DEG = 90;

  // UI microscopy mode preset (does NOT change renderer math in Wave 1)
  var mode = "TEM";

  // Provider-based state (lazy-ready)
  var provider = null;
  var providerMeta = null;
  var title = "";

  // last chosen source
  var active_source = "smiles"; // 'smiles' | 'file'

  function updateActiveSourceIndicator() {
    // Placeholder for active input indicator (SMILES vs file).
  }

  // SMILES error de-dup (avoid alert spam for same input)
  var lastSmilesErrorToken = "";

  // SMILES atom identity cache (fixes cases where SMILES provider returns wrong/empty Z/sym).
  // Filled on successful SMILES build; applied before rendering + element overrides table.
  var _smilesAtomZCache = { token: "", Zs: null, syms: null };

  // currently loaded file (if any)
  var file_state = {
    kind: null, // 'cif' | 'csv' | 'xyz' | 'json' | 'poscar' | null
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

  // Centralized SMILES rebuild helper (used by input, Enter, and Samples dropdown).
  function scheduleSmilesRebuild(opts) {
    active_source = "smiles";
    setSmilesErrorVisible(false);
    lastSmilesErrorToken = "";
    var explicit = opts && opts.explicit === true;
    rebuildAndRender({ explicit: explicit });
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

  function collectSettingsFromUI() {
    var cfg = {};
    cfg.mode = mode;
    cfg.activeSource = active_source;

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
      active_source = "smiles";
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
      active_source = "smiles";
      updateActiveSourceIndicator();
      setSmilesErrorVisible(false);
      scheduleSmilesRebuild({ explicit: true });
      setSamplesStatus(null);
      if (selSamples) selSamples.value = "__placeholder__";
      saveSettingsDebounced();
      return;
    }

    if (
      sample.kind === "file" &&
      sample.path &&
      window.fetch &&
      typeof handleLocalFile === "function"
    ) {
      setSamplesStatus("sampleLoading");
      try {
        var resp = await fetch(sample.path);
        if (!resp.ok) {
          throw new Error("HTTP " + resp.status);
        }
        var text = await resp.text();
        var fileName =
          sample.name || sample.path.split("/").pop() || "sample.dat";
        var blob = new Blob([text], { type: "text/plain" });
        var file = new File([blob], fileName, { type: "text/plain" });
        await handleLocalFile(file);
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

  function isPeriodicFileKind(kind, name) {
    if (kind === "cif" || kind === "json" || kind === "poscar") return true;
    var low = String(name || "").toLowerCase();
    if (low === "poscar" || low === "contcar") return true;
    return false;
  }

  // ---- build + render ----
  async function buildCurrentSystem(cfg) {
    cfg = cfg || {};
    var explicit = !!cfg.explicit;

    var spec = "";
    var kind = null;
    var isSmiles = false;
    var smilesToken = "";

    if (active_source === "file" && file_state.url) {
      // when loading files, hide SMILES error (not relevant)
      setSmilesErrorVisible(false);

      // extension is encoded in #fragment so phases can detect .cif/.csv
      spec = file_state.url + "#" + file_state.name;
      kind = file_state.kind;

      if (isPeriodicFileKind(kind, file_state.name)) {
        var tiles = parseTilesSize();
        if (tiles) spec = applyPeriodicSizeToSpec(spec, tiles);
      }
    } else {
      active_source = "smiles";
      smilesToken = tbSmiles ? trimStr(tbSmiles.value) : "";
      spec = smilesToken || "O";
      kind = "smiles";
      isSmiles = true;
    }

    // SMILES: handle invalid input without breaking the app or overwriting the last good image
    if (isSmiles) {
      var msg = getSmilesInvalidMsg();

      function failSmiles() {
        setSmilesErrorVisible(true, msg);
        if (explicit) {
          maybeAlertSmilesError(msg, smilesToken, true);
        }
        return false; // keep previous provider/title/canvas
      }

      var newProvider = null;
      try {
        var _smOk = true;
        if (smilesToken) {
          _smOk = await rdkitValidateSmiles(smilesToken);
          if (!_smOk) return failSmiles();
        }

        newProvider = await build_system_from_input(spec);
      } catch (e) {
        return failSmiles();
      }

      // Some parsers may fail silently and produce empty atoms.
      // Treat empty atoms as invalid ONLY when the user actually typed a SMILES string.
      if (smilesToken) {
        try {
          var testViewState = make_view_state({
            img_size: view.img_size,
            angstroms_per_pixel: view.angstroms_per_pixel,
            pan_px: view.pan_px,
            rotZ_rad: (view.rot_deg * Math.PI) / 180.0,
            center_A: view.center_A || [0, 0, 0],
            center_mode: "bbox",
          });
          var provView =
            newProvider && typeof newProvider.getView === "function"
              ? newProvider.getView(testViewState, { needBonds: false })
              : null;
          var atomsTest =
            provView && provView.atomsView ? provView.atomsView : [];
          if (!atomsTest || atomsTest.length === 0) return failSmiles();
        } catch (e2) {
          return failSmiles();
        }
      }

      provider = newProvider;

      // Cache per-atom identity from RDKit for SMILES to avoid cases where provider atoms lose Z/sym.
      try {
        var _smilesAtomInfo = await rdkitGetAtomZs(spec);
        _smilesAtomZCache.token = spec;
        _smilesAtomZCache.Zs =
          _smilesAtomInfo && Array.isArray(_smilesAtomInfo.Zs)
            ? _smilesAtomInfo.Zs
            : null;
        _smilesAtomZCache.syms =
          _smilesAtomInfo && Array.isArray(_smilesAtomInfo.syms)
            ? _smilesAtomInfo.syms
            : null;
      } catch (e3) {
        _smilesAtomZCache.token = spec;
        _smilesAtomZCache.Zs = null;
        _smilesAtomZCache.syms = null;
      }

      // success => hide error and reset de-dup token
      setSmilesErrorVisible(false);
      lastSmilesErrorToken = "";
    } else {
      provider = await build_system_from_input(spec);
    }

    providerMeta =
      provider && typeof provider.getMeta === "function"
        ? provider.getMeta()
        : null;
    if (providerMeta && Array.isArray(providerMeta.center_A)) {
      // Keep a stable anchor for pan/rot and for future lazy-tiling selection.
      view.center_A = providerMeta.center_A;
    }
    var uiTitle =
      active_source === "file" && file_state && file_state.name
        ? file_state.name
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

    if (!provider || !providerMeta) {
      clearCanvas(getBackgroundGray());
      return;
    }

    // keep viewState in sync with current render params
    view.canvas_w = w;
    view.canvas_h = h;
    view.img_size = [h, w];
    view.angstroms_per_pixel = rngZoom
      ? parseFloat(rngZoom.value)
      : view.angstroms_per_pixel;

    // build viewState early (camera affects BOTH selection and projection)
    var viewState = make_view_state({
      img_size: view.img_size,
      angstroms_per_pixel: view.angstroms_per_pixel,
      pan_px: view.pan_px,
      rotZ_rad: (view.rot_deg * Math.PI) / 180.0,
      center_A: view.center_A || [0, 0, 0],
      center_mode: "bbox",
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

    var provView = provider.getView(viewState, {
      needBonds: cbBonds ? !!cbBonds.checked : true,
    });
    var atoms = provView && provView.atomsView ? provView.atomsView : [];
    var bonds = provView && "bondsView" in provView ? provView.bondsView : null;
    var usedCamera =
      provView && provView.usedCamera ? provView.usedCamera : viewState;
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

    // SMILES: ensure atoms have correct Z (and sym) before overrides/UI gather.
    if (active_source === "smiles") {
      try {
        applySmilesAtomZsIfAvailable(atoms);
      } catch (e0) {}
    }

    // Update element overrides UI (after we have atoms for this view)
    try {
      var elems = gatherElementsFromAtoms(atoms);
      var key = elems.join("|");
      if (key !== _elemListKey) {
        _elemListKey = key;
        rebuildElementOverridesTable(elems);
      }
    } catch (e) {
      // silent
    }

    // focal_z from atoms
    var focus_norm = rngFocus ? parseFloat(rngFocus.value) : 0.5;
    if (!Number.isFinite(focus_norm)) focus_norm = 0.5;

    // focal_z from atoms (respect model tilt X/Y if used)
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
    var focal_z = zmin + focus_norm * (zmax - zmin);

    var opts = {
      bonds: bonds,
      img_size: [h, w],
      angstroms_per_pixel: viewState.angstroms_per_pixel,
      blur_sigma: rngBlur ? parseFloat(rngBlur.value) : 1.0,
      background_gray: getBackgroundGray(),
      invert: cbInvert ? !!cbInvert.checked : false,
      noise_stddev:
        cbNoise && cbNoise.checked && rngNoise
          ? parseFloat(rngNoise.value)
          : 0.0,
      contrast: rngContrast ? parseFloat(rngContrast.value) : 1.0,
      compose_mode: "sum",
      draw_bonds_flag: cbBonds ? !!cbBonds.checked : true,
      camera: usedCamera,
      bond_wave_width_px: rngBwidth ? parseFloat(rngBwidth.value) : 6,
      bond_wave_amplitude: rngBamp ? parseFloat(rngBamp.value) : 0.4,
      low_clip:
        el("tb-clip-lo") && el("tb-clip-lo").value !== ""
          ? parseFloat(el("tb-clip-lo").value)
          : null,
      high_clip:
        el("tb-clip-hi") && el("tb-clip-hi").value !== ""
          ? parseFloat(el("tb-clip-hi").value)
          : null,
      focal_z: focal_z,
      dof_strength: mode === "TEM" && rngDof ? parseFloat(rngDof.value) : 0.0,
      hide_front: cbHide ? !!cbHide.checked : false,
      show_scale_bar: cbScale ? !!cbScale.checked : false,
      scale_bar_corner: "bl",
      scale_bar_margin_px: 12,
      canvasCtx: ctx,

      // UI-only (Wave 1): renderer may ignore these for now
      mode: mode,
      scanlines: cbScanlines ? !!cbScanlines.checked : false,
      tip_sharpness: rngTip ? parseFloat(rngTip.value) : 0.5,

      element_overrides: elementOverrides,
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
  async function handleLocalFile(file) {
    if (!file) return;
    var kind = fileKindFromName(file.name);
    // Allow extensionless POSCAR/CONTCAR and content-sniffing in phases.js (blob urls).
    // Keep a light safety net: ignore huge binary files.
    if (!kind && (file.size | 0) > 20 * 1024 * 1024) return; // 20 MB

    revokeFileUrl();

    file_state.kind = kind;
    file_state.name = file.name;
    file_state.url = URL.createObjectURL(file);
    file_state.isBlob = true;

    active_source = "file";
    await rebuildAndRender();
  }

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
    if (
      active_source === "file" &&
      isPeriodicFileKind(file_state.kind, file_state.name)
    )
      rebuildAndRender();
  }
  [tbNx, tbNy, tbNz].forEach(function (t) {
    if (t) t.onchange = onTilesChanged;
  });

  if (tbSmiles) {
    tbSmiles.value = trimStr(smiles_text);

    // Default SMILES on first load: keep it simple.
    if (!trimStr(tbSmiles.value)) tbSmiles.value = "O";

    // Switching from file -> SMILES must be immediate.
    // Rebuild on any change (typing, paste, samples).
    tbSmiles.oninput = function () {
      scheduleSmilesRebuild({ explicit: false });
    };
    tbSmiles.onchange = function () {
      scheduleSmilesRebuild({ explicit: false });
    };

    // Make plain Enter apply.
    // - For <textarea>: Enter applies; Shift+Enter inserts newline.
    // - For <input>: Enter applies.
    tbSmiles.addEventListener("keydown", function (ev) {
      if (ev.key !== "Enter") return;

      var isTextarea =
        String(tbSmiles.tagName || "").toUpperCase() === "TEXTAREA";
      if (isTextarea && ev.shiftKey) return; // allow newline

      ev.preventDefault();
      scheduleSmilesRebuild({ explicit: true });
    });
  }

  if (fileInput) {
    fileInput.onchange = function () {
      var files = fileInput.files;
      if (files && files.length) handleLocalFile(files[0]);
      fileInput.value = "";
    };
  }

  // Global drag&drop (drop anywhere)
  document.addEventListener("dragover", function (ev) {
    ev.preventDefault();
  });
  document.addEventListener("drop", function (ev) {
    ev.preventDefault();
    var dt = ev.dataTransfer;
    if (!dt || !dt.files || !dt.files.length) return;
    handleLocalFile(dt.files[0]);
  });

  if (btnExport) {
    btnExport.onclick = function () {
      var a = document.createElement("a");
      a.download = "tem.png";
      a.href = cvs.toDataURL("image/png");
      a.click();
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

  // View controls: pan (px) + rotation (deg) — rerender only (NO rebuild/parsing)
  function updateViewAndRender(fn) {
    if (typeof fn === "function") fn();
    scheduleRender();
  }

  if (btnViewUp)
    btnViewUp.onclick = function () {
      updateViewAndRender(function () {
        view.pan_px[1] -= PAN_STEP_PX;
      });
    };
  if (btnViewDown)
    btnViewDown.onclick = function () {
      updateViewAndRender(function () {
        view.pan_px[1] += PAN_STEP_PX;
      });
    };
  if (btnViewLeft)
    btnViewLeft.onclick = function () {
      updateViewAndRender(function () {
        view.pan_px[0] -= PAN_STEP_PX;
      });
    };
  if (btnViewRight)
    btnViewRight.onclick = function () {
      updateViewAndRender(function () {
        view.pan_px[0] += PAN_STEP_PX;
      });
    };
  if (btnViewRotL)
    btnViewRotL.onclick = function () {
      updateViewAndRender(function () {
        view.rot_deg -= ROT_STEP_DEG;
      });
    };
  if (btnViewRotR)
    btnViewRotR.onclick = function () {
      updateViewAndRender(function () {
        view.rot_deg += ROT_STEP_DEG;
      });
    };
  if (btnViewReset)
    btnViewReset.onclick = function () {
      updateViewAndRender(function () {
        view.pan_px = [0, 0];
        view.rot_deg = 0;
      });
    };

  // Model tilt: rotate structure itself by +90° around X or Y
  if (btnModelX90)
    btnModelX90.onclick = function () {
      updateViewAndRender(function () {
        view.tilt_x_deg = (view.tilt_x_deg + TILT_STEP_DEG) % 360;
      });
    };
  if (btnModelY90)
    btnModelY90.onclick = function () {
      updateViewAndRender(function () {
        view.tilt_y_deg = (view.tilt_y_deg + TILT_STEP_DEG) % 360;
      });
    };

  // refresh dynamic labels on language switch
  try {
    document.addEventListener("i18n-updated", function () {
      updateGifUI();
      updateModeDescText();
      updateBondsLabelText();
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
      clearSettingsStorage();
      // Reset to clean TEM defaults; HTML already encodes default slider values.
      selMode.value = "TEM";
      applyModePreset("TEM", { save: false, rerender: false });
      // Re-render with defaults and store them again.
      scheduleRender();
    });
  }

  // Default active source is SMILES; in restored state this is already true
  // for cases where we can actually rebuild content.
  active_source = "smiles";
  getBackgroundGray();
  updateGifUI();

  // If we have SMILES text after restore, rebuild from it; otherwise just draw empty view
  if (tbSmiles && tbSmiles.value) {
    await rebuildAndRender();
  } else {
    scheduleRender();
  }
}
