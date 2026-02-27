// js/chem/rdkit_wrap.js
// Безпечна обгортка для RDKit_minimal.js, що працює навіть при затримці ініціалізації WASM

let _RDKit = null; // глобальний екземпляр RDKitModule

// --- внутрішня функція отримання модуля ---
function req() {
    if (!_RDKit) throw new Error("RDKitModule ще не ініціалізований");
    return _RDKit;
}

// --- головний фасад RDKit ---
export const RDKit = {
    setModule(mod) {
        _RDKit = mod;
    },

    get_mol(sm) {
        const mod = req();
        if (typeof mod.get_mol === "function") return mod.get_mol(sm);
        if (mod.Mol) return new mod.Mol(sm);
        throw new Error("RDKit: get_mol недоступна у цій збірці");
    },

    ptable_getAtomicNumber(sym) {
        const mod = req();
        return mod.getAtomicNumber ? mod.getAtomicNumber(sym) : 0;
    },

    ptable_getRcovalent(Z) {
        const mod = req();
        return mod.getRcovalent ? mod.getRcovalent(Z) : 0.77;
    },

    // координати через MolBlock (fallback)
    get_coordinates(mol) {
        const mb = mol?.get_molblock?.() ?? mol?.MolToMolBlock?.();
        if (!mb) return [];
        const lines = mb.split("\n");

        // Counts line (рядок 4, index=3) формату V2000
        // Перші 3 колонки — кількість атомів, наступні 3 — кількість зв'язків
        const counts = lines[3] ?? "";
        const nat = parseInt(counts.slice(0, 3), 10) || 0;

        const atoms = [];
        for (let i = 0; i < nat; i++) {
            const ln = lines[4 + i] ?? "";
            const x = parseFloat(ln.slice(0, 10));
            const y = parseFloat(ln.slice(10, 20));
            const z = parseFloat(ln.slice(20, 30));
            if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
                atoms.push({ x, y, z });
            }
        }
        return atoms;
    },

    get_atoms_bonds(mol) {
        const mod = req();
        const mb = mol?.get_molblock?.() ?? mol?.MolToMolBlock?.();
        if (!mb) return { atoms: [], bonds: [] };
        const lines = mb.split("\n");

        const counts = lines[3] ?? "";
        const nat = parseInt(counts.slice(0, 3), 10) || 0;
        const nb = parseInt(counts.slice(3, 6), 10) || 0;

        const atoms = [];
        for (let i = 0; i < nat; i++) {
            const ln = lines[4 + i] ?? "";
            const x = parseFloat(ln.slice(0, 10));
            const y = parseFloat(ln.slice(10, 20));
            const z = parseFloat(ln.slice(20, 30));
            const sym = ln.slice(31, 34).trim(); // колонки 31–34 — символ елемента
            const Z = (mod.getAtomicNumber ? mod.getAtomicNumber(sym) : 0) || 6; // фолбек на C
            atoms.push({ Z, x, y, z });
        }

        const bonds = [];
        // блок зв'язків починається після атомних рядків:
        const bondStart = 4 + nat;
        for (let k = 0; k < nb; k++) {
            const ln = lines[bondStart + k] ?? "";
            const i = (parseInt(ln.slice(0, 3), 10) || 0) - 1; // molfile індексує з 1
            const j = (parseInt(ln.slice(3, 6), 10) || 0) - 1;
            const order = parseInt(ln.slice(6, 9), 10) || 1;  // 1,2,3 (ароматика в molfile теж може бути 4)
            if (i >= 0 && j >= 0) bonds.push([i, j, order]);
        }

        return { atoms, bonds };
    }


};

// --- bootstrap ініціалізація ---
const rdkitReadyPromise = (() => {
    if (window.RDKitReady && typeof window.RDKitReady.then === "function") {
        return window.RDKitReady;
    }
    if (typeof window.initRDKitModule === "function") {
        const p = window.initRDKitModule({
            locateFile: (path) => "./rdkit/RDKit_minimal.wasm"
        });
        window.RDKitReady = p;
        return p;
    }
    throw new Error("RDKit_minimal.js не підключений або підключений занадто пізно.");
})();

// --- чекаємо готовності RDKit ---
rdkitReadyPromise.then(mod => {
    _RDKit = mod; // тепер модуль гарантовано готовий

    // якщо є клас Mol — додаємо допоміжні методи
    if (mod.Mol && mod.Mol.prototype) {
        const M = mod.Mol.prototype;

        M.add_hs = function () { try { this.addHs?.(); } catch { } };
        M.embed_molecule = function () { try { return this.embedMolecule?.() ?? false; } catch { return false; } };
        M.uff_optimize_molecule = function () { try { this.UFFOptimizeMolecule?.(); } catch { } };
        M.generate_aligned_coords = function () { try { this.generate_alignedCoords?.(); } catch { } };

        M.get_num_atoms = function () { try { return this.getNumAtoms?.() ?? 0; } catch { return 0; } };
        M.get_atom_atomicnum = function (i) { try { return this.getAtomWithIdx?.(i)?.getAtomicNum?.() ?? 6; } catch { return 6; } };

        M.get_num_bonds = function () { try { return this.getNumBonds?.() ?? 0; } catch { return 0; } };
        M.get_bond_begin = function (k) { try { return this.getBondWithIdx?.(k)?.getBeginAtomIdx?.() ?? 0; } catch { return 0; } };
        M.get_bond_end = function (k) { try { return this.getBondWithIdx?.(k)?.getEndAtomIdx?.() ?? 0; } catch { return 0; } };
        M.get_bond_order = function (k) {
            try {
                const b = this.getBondWithIdx?.(k);
                if (!b) return 1;
                if (b.getIsAromatic?.()) return 4;
                return b.getBondTypeAsDouble?.() ?? 1;
            } catch { return 1; }
        };
    }

    console.log("✅ RDKit WASM ініціалізовано");
}).catch(e => {
    console.error("❌ Помилка ініціалізації RDKit:", e);
});
