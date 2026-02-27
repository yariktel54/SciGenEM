// js/io/smiles_io.js
// create_molecule_from_smiles(smiles), get_atoms_with_coords(mol)
import { RDKit } from '../chem/rdkit_wrap.js';
import { is_dev_mode } from '../system/system_contract.js';

// Розкласти незв’язані компоненти (якщо RDKit не дав 2D-розкладку)
function separate_fragments_inplace(atoms, bonds, gap = 8.0) {
    if (!atoms.length) return;

    // Будуємо граф за bonds
    const g = new Map();
    for (let i = 0; i < atoms.length; i++) g.set(i, []);
    for (const [i, j] of bonds) {
        if (i >= 0 && j >= 0 && i < atoms.length && j < atoms.length) {
            g.get(i).push(j);
            g.get(j).push(i);
        }
    }

    // Зв’язні компоненти (DFS)
    const comps = [];
    const seen = new Array(atoms.length).fill(false);
    const stack = [];
    for (let s = 0; s < atoms.length; s++) {
        if (seen[s]) continue;
        stack.length = 0;
        stack.push(s);
        seen[s] = true;
        const comp = [];
        while (stack.length) {
            const v = stack.pop();
            comp.push(v);
            for (const u of g.get(v)) {
                if (!seen[u]) { seen[u] = true; stack.push(u); }
            }
        }
        comps.push(comp);
    }
    if (comps.length <= 1) return; // нема що розводити

    // Центруємо кожну компоненту навколо (0,0) і розносимо по X з кроком gap
    let shiftX = 0;
    for (const comp of comps) {
        const cx = comp.reduce((s, i) => s + atoms[i].x, 0) / comp.length;
        const cy = comp.reduce((s, i) => s + atoms[i].y, 0) / comp.length;
        for (const i of comp) {
            atoms[i].x = (atoms[i].x - cx) + shiftX;
            atoms[i].y = (atoms[i].y - cy);
        }
        shiftX += gap;
    }
}


// --- TEM helpers: non-fatal warnings ---
function _tem_push_warning(mol, msg) {
    try {
        if (!mol) return;
        if (!Array.isArray(mol.__tem_warnings)) mol.__tem_warnings = [];
        mol.__tem_warnings.push(String(msg));
    } catch (_) { /* ignore */ }
}

function _tem_is_degenerate_coords(atoms, eps = 1e-6) {
    if (!Array.isArray(atoms) || atoms.length < 2) return true;
    let minx = atoms[0].x, maxx = atoms[0].x;
    let miny = atoms[0].y, maxy = atoms[0].y;
    for (let i = 1; i < atoms.length; i++) {
        const a = atoms[i];
        if (!Number.isFinite(a.x) || !Number.isFinite(a.y)) return true;
        if (a.x < minx) minx = a.x; if (a.x > maxx) maxx = a.x;
        if (a.y < miny) miny = a.y; if (a.y > maxy) maxy = a.y;
    }
    return (maxx - minx) < eps && (maxy - miny) < eps;
}

function _tem_apply_fallback_2d_layout(atoms, r = 3.0) {
    const n = atoms.length;
    if (!n) return;
    for (let i = 0; i < n; i++) {
        const t = (2 * Math.PI * i) / n;
        atoms[i].x = Math.cos(t) * r;
        atoms[i].y = Math.sin(t) * r;
        atoms[i].z = 0;
    }
}

export async function create_molecule_from_smiles(smiles) {
    const DEV = is_dev_mode();
    const dlog = (...a) => { if (DEV) console.debug(...a); };
    const dwarn = (...a) => { if (DEV) dwarn(...a); };
    dlog("=== create_molecule_from_smiles start ===");
    dlog("SMILES:", smiles);

    // 1. Гарантуємо готовність RDKit
    if (window.RDKitReady && typeof window.RDKitReady.then === "function") {
        const mod = await window.RDKitReady;
        if (!window.RDKit) {
            const wrap = await import("../chem/rdkit_wrap.js");
            window.RDKit = wrap.RDKit;
            dlog("RDKit wrap імпортовано у window.");
        }
        if (typeof window.RDKit.setModule === "function") {
            window.RDKit.setModule(mod);
            dlog("RDKit.setModule виконано.");
        }
    } else {
        dwarn("❌ RDKitReady не знайдено або ще не завантажений.");
        throw new Error("RDKit не готовий");
    }

    // 2. Створюємо молекулу
    const mol = window.RDKit.get_mol(smiles);
    try { if (!Array.isArray(mol.__tem_warnings)) mol.__tem_warnings = []; } catch (_) {}
    if (!mol) throw new Error("Не вдалося зчитати SMILES: " + smiles);
    dlog("Mol створено:", mol);

    // 3. Кількість атомів/зв’язків одразу після створення
    const nat_init = mol.get_num_atoms?.() ?? -1;
    const nb_init = mol.get_num_bonds?.() ?? -1;
    dlog(`  → атомів = ${nat_init}, зв’язків = ${nb_init}`);

    // 4. Додаємо H
    try { mol.add_hs?.(); dlog("  + add_hs виконано"); } catch { dwarn("  add_hs недоступна"); }

    // 5. Спроба створити 3D
    let ok = false;
    try {
        ok = mol.embed_molecule?.() || mol.EmbedMolecule?.() || false;
        dlog("  embed_molecule:", ok);
        if (ok) {
            mol.uff_optimize_molecule?.();
            mol.UFFOptimizeMolecule?.();
            dlog("  UFF optimize виконано");
        }
    } catch (e) {
        dwarn("  embed/optimize виняток", e);
        ok = false;
    }

    // 6. Якщо ні — 2D (RDKit_minimal інколи вимагає аргумент)
    if (!ok) {
        let did2d = false;
        // Prefer canonical 2D APIs, if present
        try {
            if (typeof mol.compute2DCoords === 'function') { mol.compute2DCoords(); did2d = true; }
            else if (typeof mol.compute_2d_coords === 'function') { mol.compute_2d_coords(); did2d = true; }
        } catch (e) {
            _tem_push_warning(mol, 'SMILES: compute2DCoords failed');
            dwarn('  compute2DCoords виняток', e);
        }

        // Fallback: generate_aligned_coords may require a reference mol argument
        if (!did2d) {
            try {
                if (typeof mol.generate_aligned_coords === 'function') { mol.generate_aligned_coords(mol); did2d = true; }
                else if (typeof mol.generate_alignedCoords === 'function') { mol.generate_alignedCoords(mol); did2d = true; }
            } catch (e) {
                _tem_push_warning(mol, 'SMILES: generate_aligned_coords failed');
                dwarn('  generate_aligned_coords виняток', e);
                // Last chance: retry without args
                try {
                    mol.generate_aligned_coords?.();
                    mol.generate_alignedCoords?.();
                    did2d = true;
                } catch (_) {}
            }
        }
    }
// 7. Отримуємо координати
    const coords = window.RDKit.get_coordinates?.(mol) ?? [];
    dlog("  get_coordinates:", coords.length, coords);

    // 8. Якщо порожньо — створюємо випадкові
    if (!coords.length) {
        const nat = mol.get_num_atoms?.() ?? 0;
        const fake = [];
        for (let i = 0; i < nat; i++) {
            fake.push({ x: Math.random() * 10 - 5, y: Math.random() * 10 - 5, z: 0 });
        }
        mol.get_coordinates = () => fake;
        _tem_push_warning(mol, 'SMILES: coords missing; used random fallback');
        dwarn("⚠️ координати відсутні — створено випадкові позиції:", fake);
    }

    // 9. Підсумок
    dlog("Mol atoms:", mol.get_num_atoms?.(), "bonds:", mol.get_num_bonds?.());
    dlog("=== create_molecule_from_smiles end ===");

    return mol;
}


export function get_atoms_with_coords(mol, use_2d_if_fail = true) {
    // Спочатку намагаємось витягнути все з MolBlock
    if (typeof window.RDKit?.get_atoms_bonds === "function") {
        const { atoms, bonds } = window.RDKit.get_atoms_bonds(mol);
        if (atoms.length) {
            // Якщо координати вироджені (усі ~0), застосувати простий 2D-layout.
            if (_tem_is_degenerate_coords(atoms)) {
                _tem_push_warning(mol, 'SMILES: degenerate coords; applied fallback 2D layout');
                _tem_apply_fallback_2d_layout(atoms, 3.0);
            }
            // Якщо зв’язків мало щодо атомів — ймовірно кілька фрагментів → розвести
            if (bonds.length < atoms.length - 1) {
                separate_fragments_inplace(atoms, bonds, 8.0); // 8 Å між компонентами
            }
            // Тепер фінальне центровання всієї системи
            const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length;
            const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length;
            for (const a of atoms) { a.x -= cx; a.y -= cy; }
            return [atoms, bonds];
        }
    }


    // Інакше — старий шлях (координати без Z/бonds → фолбек)
    const conf = window.RDKit.get_coordinates?.(mol) ?? [];
    if (!conf.length) {
        if (use_2d_if_fail) {
            try { mol.generate_aligned_coords?.(); mol.generate_alignedCoords?.(); } catch { }
            return get_atoms_with_coords(mol, false);
        }
        throw new Error("❌ Молекула не має координат.");
    }

    const natoms = conf.length;
    const atoms = [];
    for (let i = 0; i < natoms; i++) {
        // фолбек: якщо не змогли дістати Z — C
        const Z = mol.get_atom_atomicnum?.(i) ?? 6;
        const { x, y, z } = conf[i];
        atoms.push({ Z, x, y, z });
    }
    const bonds = []; // без MolBlock — зв’язків не дістанемо; renderer їх здогадається по відстані
    // Центрування
    const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length;
    const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length;
    for (const a of atoms) { a.x -= cx; a.y -= cy; }

    return [atoms, bonds];
}


