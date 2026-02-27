// js/app/phases.js
// System builder: detects input format (SMILES / CIF / CSV / XYZ / JSON / MOL / PDB) and returns atoms + bonds.
// NOTE: "liquid/solid" modes were removed. We always build what the last input represents.
//
// Safe-refactor wave A:
//  - phases.js is a thin dispatcher
//  - CIF tiling math moved to lattice.js
//  - per-format building moved to builders/

import { build_from_smiles } from '../builders/build_smiles.js';
import { build_from_cif } from '../builders/build_cif.js'; // kept for compatibility; CIF provider bypasses eager tiling
import { build_from_csv } from '../builders/build_csv.js';
import { build_from_xyz } from '../builders/build_xyz.js';
import { build_from_json } from '../builders/build_json.js';
import { build_from_poscar } from '../builders/build_poscar.js';
import { build_from_mol } from '../builders/build_mol.js';
import { build_from_pdb } from '../builders/build_pdb.js';

import { is_dev_mode, normalize_build_result, validate_system } from '../system/system_contract.js';

import { load_cif } from '../io/cif_io.js';
import { load_json_system } from '../io/json_io.js';
import { StaticProvider, PeriodicProvider } from '../system/system_provider.js';
import { cell_vectors, frac_to_cart } from '../system/lattice.js';

export function parse_system_string(s) {
    s = (s || '').trim();
    const low = s.toLowerCase();

    // ---------- helpers ----------
    const parse_size = (val) => {
        val = (val || '').trim().toLowerCase();
        if (val === 'auto') return 'auto';
        const parts = val.split(/[x×]/);
        if (parts.length !== 3) return [1, 1, 1];
        const nx = parseInt(parts[0], 10);
        const ny = parseInt(parts[1], 10);
        const nz = parseInt(parts[2], 10);
        if ([nx, ny, nz].some(Number.isNaN)) return [1, 1, 1];
        return [nx, ny, nz];
    };

    const parse_kv = (body) => {
        const items = body.replace(/,/g, ';').split(';').map(x => x.trim()).filter(Boolean);
        const kv = {};
        for (const it of items) {
            if (it.includes('=')) {
                const [k, v] = it.split('=', 2);
                kv[k.trim().toLowerCase()] = (v ?? '').trim();
            }
        }
        return kv;
    };

    const looks_like_csv_text = (txt) => {
        // Strict CSV uses section header "ATOMS" and optional "BONDS".
        return /(^|\n)\s*ATOMS\s*(\n|$)/i.test(txt);
    };

    const looks_like_cif_text = (txt) => {
        // CIF usually contains data_ and/or atom_site/cell fields.
        if (/(^|\n)\s*data_/i.test(txt)) return true;
        if (/(^|\n)\s*_atom_site_/i.test(txt)) return true;
        if (/(^|\n)\s*_cell_length_/i.test(txt)) return true;
        return false;
    };

    const looks_like_xyz_text = (txt) => {
        // XYZ: first line is integer N; must be multi-line.
        const s0 = String(txt || '');
        if (!s0.includes('\n')) return false;
        const lines = s0.split(/\r?\n/);
        if (!lines.length) return false;
        const n = parseInt(String(lines[0] || '').trim(), 10);
        if (!Number.isFinite(n) || n <= 0) return false;
        // be conservative: require at least N atom lines + comment line
        return lines.length >= (n + 2);
    };

    const looks_like_poscar_text = (txt) => {
        // POSCAR/CONTCAR (VASP5) heuristic: scale on line2, then 3 lattice vectors.
        const s0 = String(txt || '');
        if (!s0.includes('\n')) return false;
        const lines = s0.split(/\r?\n/).map(x => String(x || '').trim()).filter(x => x !== '');
        if (lines.length < 8) return false;
        const scale = parseFloat(lines[1]);
        if (!Number.isFinite(scale) || scale === 0) return false;
        for (let i = 2; i <= 4; i++) {
            const p = lines[i].split(/\s+/).filter(Boolean);
            if (p.length < 3) return false;
            const v = [parseFloat(p[0]), parseFloat(p[1]), parseFloat(p[2])];
            if (!v.every(Number.isFinite)) return false;
        }
        const symT = lines[5].split(/\s+/).filter(Boolean);
        const cntT = lines[6].split(/\s+/).filter(Boolean);
        if (!symT.length || !cntT.length) return false;
        if (!symT.every(t => /^[A-Za-z]{1,2}$/.test(t))) return false;
        if (!cntT.every(t => /^-?\d+$/.test(t))) return false;
        return true;
    };

    const looks_like_json_text = (txt) => {
        // JSON: starts with { or [ (after trimming / BOM). Avoid mis-detecting SMILES like "[Ni+2]".
        let t = String(txt || '');
        if (!t) return false;
        t = t.replace(/^\uFEFF/, '').trim();
        if (!t) return false;
        const c0 = t[0];
        if (!(c0 === '{' || c0 === '[')) return false;
        // Require a likely key to avoid collisions with bracketed SMILES.
        return /"atoms"\s*:/i.test(t) || /"bonds"\s*:/i.test(t) || /"cell"\s*:/i.test(t) || /"kind"\s*:/i.test(t) || /"sites"\s*:/i.test(t) || /"lattice"\s*:/i.test(t);
    };

    const looks_like_mol_text = (txt) => {
        // MOL: contains counts line with V2000/V3000 and ends with 'M  END'.
        const s0 = String(txt || '');
        if (!s0.includes('\n')) return false;
        if (/(^|\n)M\s{2}END(\s|$)/.test(s0)) return true;
        const lines = s0.split(/\r?\n/);
        if (lines.length >= 4) {
            const counts = String(lines[3] || '');
            if (/V2000|V3000/i.test(counts)) return true;
        }
        return false;
    };

    const looks_like_pdb_text = (txt) => {
        const s0 = String(txt || "");
        if (!s0.includes("\n")) return false;
        // PDB: fixed-record format. Look for common record types.
        if (/(^|\n)ATOM  /i.test(s0)) return true;
        if (/(^|\n)HETATM/i.test(s0)) return true;
        if (/(^|\n)HEADER/i.test(s0)) return true;
        return false;
    };

    // -------- format registry: add new formats by appending handlers --------
    // Each handler: { name, detect(str, low, helpers) -> spec|null }

    const helpers = { parse_size, parse_kv, looks_like_csv_text, looks_like_cif_text, looks_like_xyz_text, looks_like_json_text, looks_like_poscar_text, looks_like_mol_text, looks_like_pdb_text };

    const detect_cif = (str, lowStr, h) => {
        // CIF with optional args ("...cif;size=2x2x1")
        // Works for normal URLs/paths AND for blob urls with fragments: "blob:...#file.cif;size=...".
        if (lowStr.includes('.cif;')) {
            const [path, rest] = str.split(';', 2);
            let size = [1, 1, 1];
            const kv = h.parse_kv(rest);
            if ('size' in kv) size = h.parse_size(kv.size);
            return { mode: 'cif', path: path.trim(), size };
        }

        // Plain .cif (also matches "...#file.cif")
        if (lowStr.endsWith('.cif') || /\.cif([?#].*)?$/i.test(str)) {
            return { mode: 'cif', path: str, size: [1, 1, 1] };
        }

        // Raw pasted CIF (robustness)
        if (h.looks_like_cif_text(str)) return { mode: 'cif', path: str, size: [1, 1, 1] };
        return null;
    };

    const detect_csv = (str, lowStr, h) => {
        // Plain .csv (also matches "...#file.csv")
        if (lowStr.endsWith('.csv') || /\.csv([?#].*)?$/i.test(str)) {
            return { mode: 'csv', path: str };
        }

        // Raw pasted strict CSV (robustness)
        if (h.looks_like_csv_text(str)) return { mode: 'csv', path: str };
        return null;
    };

    const detect_xyz = (str, lowStr, h) => {
        // Plain .xyz (also matches "...#file.xyz")
        if (lowStr.endsWith('.xyz') || /\.xyz([?#].*)?$/i.test(str)) {
            return { mode: 'xyz', path: str };
        }
        // Raw pasted XYZ
        if (h.looks_like_xyz_text(str)) return { mode: 'xyz', path: str };
        return null;
    };

    const detect_pdb = (str, lowStr, h) => {
        // Plain .pdb/.ent (also matches ...#file.pdb)
        if (lowStr.endsWith(".pdb") || /\.pdb([?#].*)?$/i.test(str)) {
            return { mode: "pdb", path: str };
        }
        if (lowStr.endsWith(".ent") || /\.ent([?#].*)?$/i.test(str)) {
            return { mode: "pdb", path: str };
        }
        // Raw pasted PDB
        if (h.looks_like_pdb_text(str)) return { mode: "pdb", path: str };
        return null;
    };

    const detect_mol = (str, lowStr, h) => {
        // Plain .mol (also matches ...#file.mol)
        if (lowStr.endsWith('.mol') || /\.mol([?#].*)?$/i.test(str)) {
            return { mode: 'mol', path: str };
        }
        // Raw pasted molblock
        if (h.looks_like_mol_text(str)) return { mode: 'mol', path: str };
        return null;
    };

    const detect_poscar = (str, lowStr, h) => {
        // POSCAR/CONTCAR with optional args ("...poscar;size=2x2x1")
        if (lowStr.includes('.poscar;') || lowStr.includes('.contcar;') || /#(poscar|contcar);/i.test(str)) {
            const [path, rest] = str.split(';', 2);
            let size = [1, 1, 1];
            const kv = h.parse_kv(rest);
            if ('size' in kv) size = h.parse_size(kv.size);
            return { mode: 'poscar', path: path.trim(), size };
        }

        // Extension-based
        if (lowStr.endsWith('.poscar') || /\.poscar([?#].*)?$/i.test(str)) {
            return { mode: 'poscar', path: str, size: [1, 1, 1] };
        }
        if (lowStr.endsWith('.contcar') || /\.contcar([?#].*)?$/i.test(str)) {
            return { mode: 'poscar', path: str, size: [1, 1, 1] };
        }
        if (lowStr.endsWith('.vasp') || /\.vasp([?#].*)?$/i.test(str)) {
            return { mode: 'poscar', path: str, size: [1, 1, 1] };
        }

        // Name-only (common POSCAR/CONTCAR without extension in fragments)
        // e.g. blob:...#POSCAR
        if (/(^|[#/])(poscar|contcar)(\b|$)/i.test(str)) {
            return { mode: 'poscar', path: str, size: [1, 1, 1] };
        }

        // Raw pasted POSCAR
        if (h.looks_like_poscar_text(str)) return { mode: 'poscar', path: str, size: [1, 1, 1] };
        return null;
    };

    const detect_json = (str, lowStr, h) => {
        // JSON with optional args ("...json;size=2x2x1") for periodic tiling parity with CIF.
        if (lowStr.includes('.json;')) {
            const [path, rest] = str.split(';', 2);
            let size = [1, 1, 1];
            const kv = h.parse_kv(rest);
            if ('size' in kv) size = h.parse_size(kv.size);
            return { mode: 'json', path: path.trim(), size };
        }

        // Plain .json (also matches "...#file.json")
        if (lowStr.endsWith('.json') || /\.json([?#].*)?$/i.test(str)) {
            return { mode: 'json', path: str };
        }
        // Raw pasted JSON
        if (h.looks_like_json_text(str)) return { mode: 'json', path: str };
        return null;
    };

    const detect_smiles = (str) => {
        // Default fallback.
        return { mode: 'smiles', smiles: str };
    };

    const FORMAT_HANDLERS = [
        { name: 'cif', detect: detect_cif },
        { name: 'csv', detect: detect_csv },
        { name: 'poscar', detect: detect_poscar },
        { name: 'pdb', detect: detect_pdb },
        { name: 'mol', detect: detect_mol },
        { name: 'xyz', detect: detect_xyz },
        { name: 'json', detect: detect_json },
        { name: 'smiles', detect: (str, lowStr) => detect_smiles(str, lowStr) },
    ];

    for (const hnd of FORMAT_HANDLERS) {
        const spec = hnd.detect(s, low, helpers);
        if (spec) return spec;
    }

    // should never happen due to SMILES fallback
    return { mode: 'smiles', smiles: s };
}

export async function build_system_from_input(input_str) {
    let spec = parse_system_string(input_str);

    // Content-based fallback for extensionless file drops.
    // If spec fell back to SMILES but input is a blob/http URL, sniff text and re-dispatch.
    if (spec && spec.mode === 'smiles') {
        const s = String(input_str || '').trim();
        if (/^(https?:\/\/|blob:)/i.test(s)) {
            try {
                const r = await fetch(s);
                if (r && r.ok) {
                    const txt = await r.text();

                    const looks_like_csv = /(^|\n)\s*ATOMS\s*(\n|$)/i.test(txt);
                    const looks_like_cif = /(^|\n)\s*data_/i.test(txt) || /(^|\n)\s*_atom_site_/i.test(txt) || /(^|\n)\s*_cell_length_/i.test(txt);
                    const looks_like_xyz = (() => {
                        if (!txt.includes('\n')) return false;
                        const lines = txt.split(/\r?\n/);
                        const n = parseInt(String(lines[0] || '').trim(), 10);
                        if (!Number.isFinite(n) || n <= 0) return false;
                        return lines.length >= (n + 2);
                    })();
                    const looks_like_pdb = /(^|\n)ATOM  /i.test(txt) || /(^|\n)HETATM/i.test(txt) || /(^|\n)HEADER/i.test(txt);
                    const looks_like_mol = /(^|\n)M\s{2}END(\s|$)/.test(txt) || /V2000|V3000/i.test(txt);
                    const looks_like_json = (() => {
                        let t = String(txt || '');
                        t = t.replace(/^\uFEFF/, '').trim();
                        if (!t) return false;
                        const c0 = t[0];
                        if (!(c0 === '{' || c0 === '[')) return false;
                        return /"atoms"\s*:/i.test(t) || /"sites"\s*:/i.test(t) || /"lattice"\s*:/i.test(t) || /"cell"\s*:/i.test(t) || /"kind"\s*:/i.test(t);
                    })();
                    const looks_like_poscar = (() => {
                        if (!txt.includes('\n')) return false;
                        const lines = txt.split(/\r?\n/).map(x => String(x || '').trim()).filter(x => x !== '');
                        if (lines.length < 8) return false;
                        const scale = parseFloat(lines[1]);
                        if (!Number.isFinite(scale) || scale === 0) return false;
                        for (let i = 2; i <= 4; i++) {
                            const p = lines[i].split(/\s+/).filter(Boolean);
                            if (p.length < 3) return false;
                            const v = [parseFloat(p[0]), parseFloat(p[1]), parseFloat(p[2])];
                            if (!v.every(Number.isFinite)) return false;
                        }
                        const symT = lines[5].split(/\s+/).filter(Boolean);
                        const cntT = lines[6].split(/\s+/).filter(Boolean);
                        if (!symT.length || !cntT.length) return false;
                        if (!symT.every(t => /^[A-Za-z]{1,2}$/.test(t))) return false;
                        if (!cntT.every(t => /^-?\d+$/.test(t))) return false;
                        return true;
                    })();

                    if (looks_like_cif) spec = { mode: 'cif', path: s, size: [1, 1, 1] };
                    else if (looks_like_csv) spec = { mode: 'csv', path: s };
                    else if (looks_like_poscar) spec = { mode: 'poscar', path: s, size: [1, 1, 1] };
                    else if (looks_like_pdb) spec = { mode: 'pdb', path: s };
                    else if (looks_like_mol) spec = { mode: 'mol', path: s };
                    else if (looks_like_xyz) spec = { mode: 'xyz', path: s };
                    else if (looks_like_json) spec = { mode: 'json', path: s, size: [1, 1, 1] };
                }
            } catch (_) {
                // silent: keep SMILES fallback
            }
        }
    }

    // ------------------------------
    // CIF: true-lazy path (PeriodicProvider)
    // ------------------------------
    if (spec.mode === 'cif') {
        // spec.path may be URL/blob or raw CIF text
        const cif = await load_cif(spec.path);
        const title = (spec.path && typeof spec.path === 'string') ? (spec.path.split('#').pop().split('/').pop() || 'CIF') : 'CIF';
        const size = (spec.size === 'auto') ? 'auto' : (Array.isArray(spec.size) ? spec.size : [1, 1, 1]);

        // If cell really comes from CIF => periodic allowed.
        if (cif && cif.cell && cif.cell_from === 'cif') {
            return new PeriodicProvider({
                title,
                cell: cif.cell,
                size,
                bonds: (cif.bonds === undefined ? null : cif.bonds),
                atoms_cart: cif.atoms_cart,
                atoms_frac: cif.atoms_frac
            });
        }

        // Otherwise treat as a non-periodic molecule (no tiling).
        let atoms = [];
        if (cif && Array.isArray(cif.atoms_cart) && cif.atoms_cart.length) {
            atoms = cif.atoms_cart.map(a => ({ Z: a.Z | 0, x: +a.x || 0, y: +a.y || 0, z: +a.z || 0 }));
        } else if (cif && Array.isArray(cif.atoms_frac) && cif.atoms_frac.length) {
            // Convert frac->cart using whatever cell we have (bbox/unit fallback).
            const cell = cif.cell || [1, 1, 1, 90, 90, 90];
            const vv = cell_vectors(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
            const vx = vv[0], vy = vv[1], vz = vv[2];
            atoms = cif.atoms_frac.map(a => {
                const p = frac_to_cart([+a.fx || 0, +a.fy || 0, +a.fz || 0], vx, vy, vz);
                return { Z: a.Z | 0, x: p[0], y: p[1], z: p[2] };
            });
        }

        const bonds = (cif && ('bonds' in cif)) ? cif.bonds : null;
        return new StaticProvider({ atoms, bonds, title, meta: { mode: 'cif', cell_from: cif ? cif.cell_from : 'none' } });
    }

    // ------------------------------
    // JSON: CIF-parity path
    //  - If kind=periodic and cell exists => use PeriodicProvider (same as CIF path)
    //  - Else => StaticProvider (molecule)
    // ------------------------------
    if (spec.mode === 'json') {
        const jsys = await load_json_system(spec.path);
        const title = (spec.path && typeof spec.path === 'string') ? (spec.path.split('#').pop().split('/').pop() || 'JSON') : 'JSON';
        const size = (spec.size === 'auto') ? 'auto' : (Array.isArray(spec.size) ? spec.size : [1, 1, 1]);

        // Periodic JSON should behave like CIF: lazy tiles + view-subset.
        if (jsys && jsys.kind === 'periodic' && Array.isArray(jsys.cell) && jsys.cell.length >= 6) {
            return new PeriodicProvider({
                title,
                cell: jsys.cell,
                size,
                bonds: (jsys.bonds === undefined ? null : jsys.bonds),
                atoms_cart: jsys.atoms_cart,
                atoms_frac: jsys.atoms_frac
            });
        }

        // Otherwise treat as a non-periodic molecule.
        let atoms = [];
        if (jsys && Array.isArray(jsys.atoms_cart) && jsys.atoms_cart.length) {
            atoms = jsys.atoms_cart.map(a => ({ Z: a.Z | 0, x: +a.x || 0, y: +a.y || 0, z: +a.z || 0 }));
        } else if (jsys && Array.isArray(jsys.atoms) && jsys.atoms.length) {
            atoms = jsys.atoms.map(a => ({ Z: a.Z | 0, x: +a.x || 0, y: +a.y || 0, z: +a.z || 0 }));
        } else if (jsys && Array.isArray(jsys.atoms_frac) && jsys.atoms_frac.length) {
            const cell = jsys.cell || [1, 1, 1, 90, 90, 90];
            const vv = cell_vectors(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
            const vx = vv[0], vy = vv[1], vz = vv[2];
            atoms = jsys.atoms_frac.map(a => {
                const p = frac_to_cart([+a.fx || 0, +a.fy || 0, +a.fz || 0], vx, vy, vz);
                return { Z: a.Z | 0, x: p[0], y: p[1], z: p[2] };
            });
        }

        const bonds = (jsys && ('bonds' in jsys)) ? jsys.bonds : null;
        return new StaticProvider({ atoms, bonds, title, meta: { mode: 'json', cell_from: jsys ? jsys.cell_from : 'none', kind: jsys ? jsys.kind : 'molecule' } });
    }

    // Builder registry: add new formats here (one line) and you're done.
    const BUILDERS = {
        // cif handled above (provider lazy path)
        csv: build_from_csv,
        xyz: build_from_xyz,
        json: build_from_json,
        poscar: build_from_poscar,
        mol: build_from_mol,
        pdb: build_from_pdb,
        smiles: build_from_smiles,
        // (left here only for legacy/testing)
        cif: build_from_cif,
    };

    const build_fn = BUILDERS[spec.mode];
    if (!build_fn) throw new Error('Непідтримуваний режим');

    const warnings = [];
    let res = await build_fn(spec);

    // Merge builder-provided warnings (if any).
    try {
        const mw = res?.meta?.warnings;
        if (Array.isArray(mw)) warnings.push(...mw);
    } catch (_) {}

    const sys = normalize_build_result(res);
    validate_system(sys, `build:${spec.mode}`);

    // Lightweight warnings (no console spam). UI may choose to show them.
    if (spec.mode === 'cif') {
        if (Array.isArray(spec.size) && (spec.size[0] !== 1 || spec.size[1] !== 1 || spec.size[2] !== 1)) {
            warnings.push(`CIF tiling: ${spec.size[0]}x${spec.size[1]}x${spec.size[2]}`);
        }
        if (sys.bonds === null) {
            warnings.push('CIF: зв’язки не визначені (bonds=null)');
        } else if (Array.isArray(sys.bonds) && sys.bonds.length === 0) {
            warnings.push('CIF: зв’язки не знайдено (спробуйте збільшити size)');
        }
    }
    if (spec.mode === 'csv') {
        if (Array.isArray(sys.bonds) && sys.bonds.length === 0) {
            warnings.push('CSV: bonds порожні (BONDS секція відсутня або порожня)');
        }
    }
    if (spec.mode === 'xyz') {
        if (sys.bonds === null) warnings.push('XYZ: зв’язки не задані (bonds=null)');
    }
    if (spec.mode === 'mol') {
        warnings.push('MOL: явні зв’язки (guess заборонено)');
    }
    if (spec.mode === 'pdb') {
        if (sys.bonds === null) warnings.push('PDB: CONECT відсутній (bonds=null)');
        else if (Array.isArray(sys.bonds) && sys.bonds.length === 0) warnings.push('PDB: CONECT відсутній (bonds=[])');
        else warnings.push('PDB: CONECT (явні зв’язки)');
    }
    if (spec.mode === 'poscar') {
        warnings.push('POSCAR: зв’язки не задані (bonds=null)');
    }
    if (spec.mode === 'json') {
        if (sys.bonds === null) warnings.push('JSON: зв’язки невідомі (bonds=null)');
    }

    // Periodic via builder extras (POSCAR; also keeps door open for other formats).
    // This preserves contract semantics: sys.atoms is CART array for validation,
    // while provider uses atoms_frac when available.
    if (spec.mode === 'poscar') {
        const k = (res && res._tem_kind) ? String(res._tem_kind) : 'molecule';
        const cell = (res && Array.isArray(res._tem_cell)) ? res._tem_cell : null;
        if (k === 'periodic' && cell && cell.length >= 6) {
            const title = sys.title || 'POSCAR';
            const atoms_frac = (res && Array.isArray(res._tem_atoms_frac) && res._tem_atoms_frac.length) ? res._tem_atoms_frac : null;
            const atoms_cart = (res && Array.isArray(res._tem_atoms_cart) && res._tem_atoms_cart.length) ? res._tem_atoms_cart : sys.atoms;
            const size = (spec.size === 'auto') ? 'auto' : (Array.isArray(spec.size) ? spec.size : [1, 1, 1]);

            const prov = new PeriodicProvider({
                title,
                cell,
                size,
                bonds: sys.bonds,
                atoms_cart,
                atoms_frac
            });

            if (is_dev_mode()) {
                const nb = (sys.bonds === null) ? 'null' : String(sys.bonds.length);
                console.debug(`[TEM] build poscar(periodic): atoms=${atoms_cart.length}; bonds=${nb}`);
            }
            try {
                globalThis.TEM_LAST_BUILD = { mode: 'poscar', kind: 'periodic', atoms: atoms_cart.length, bonds: sys.bonds === null ? null : sys.bonds.length, title };
            } catch (_) {}
            return prov;
        }
    }

    // Provider return (renderer-compatible).
    const prov = new StaticProvider({
        atoms: sys.atoms,
        bonds: sys.bonds,
        title: sys.title,
        meta: {
            mode: spec.mode,
            warnings: Array.from(new Set(warnings)).slice(0, 8),
            kind: (spec.mode === 'json' && res && res._tem_kind) ? String(res._tem_kind) : undefined,
            cell: (spec.mode === 'json' && res && Array.isArray(res._tem_cell)) ? res._tem_cell : undefined
        }
    });

    // Optional DEV-only single-line debug summary.
    if (is_dev_mode()) {
        const nb = (sys.bonds === null) ? 'null' : String(sys.bonds.length);
        console.debug(`[TEM] build ${spec.mode}: atoms=${sys.atoms.length}; bonds=${nb}`);
        if (warnings.length) console.debug(`[TEM] warnings: ${warnings.join(' | ')}`);
    }

    try {
        globalThis.TEM_LAST_BUILD = { mode: spec.mode, warnings, atoms: sys.atoms.length, bonds: sys.bonds === null ? null : sys.bonds.length, title: sys.title };
    } catch (_) {}

    return prov;
}
