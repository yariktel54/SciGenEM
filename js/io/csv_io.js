// js/io/csv_io.js
// load_csv_strict(path) → (atoms, bonds) як у Python
function looksLikeRawCsvText(s) {
    s = String(s || '');
    if (!s.includes('\n')) return false;
    return /(^|\n)\s*ATOMS\s*(\n|$)/i.test(s);
}
export async function load_csv_strict(pathOrText) {
    let text = pathOrText ?? "";
    const s = String(text);
    if (!looksLikeRawCsvText(s) && /^(https?:\/\/|blob:)/i.test(s)) {
        const r = await fetch(s);
        if (!r.ok) throw new Error("CSV: не вдалося завантажити (" + r.status + ")");
        text = await r.text();
    }
    const lines0 = text.split(/\r?\n/).map(s => s.trim());
    let lines = lines0.filter((_s, i) => true);

    // обрізка порожніх країв
    while (lines.length && !lines[0]) lines.shift();
    while (lines.length && !lines.at(-1)) lines.pop();
    if (!lines.length) throw new Error("CSV порожній.");

    const findLine = (val, start = 0) => {
        for (let i = start; i < lines.length; i++) if (lines[i] === val) return i;
        return -1;
    };

    const i_atoms = findLine("ATOMS", 0);
    if (i_atoms < 0) throw new Error("Відсутній маркер секції ATOMS.");

    if (i_atoms + 1 >= lines.length) throw new Error("Очікується заголовок 'id,Z,x,y,z' після ATOMS.");
    if (lines[i_atoms + 1] !== "id,Z,x,y,z") throw new Error("Невірний заголовок ATOMS.");

    const i_bonds = findLine("BONDS", i_atoms + 2);
    let atoms_rows = lines.slice(i_atoms + 2, i_bonds >= 0 ? i_bonds : lines.length).filter(Boolean);

    const atoms = [];
    const seen = new Set();
    for (const r of atoms_rows) {
        const parts = r.split(",");
        if (parts.length !== 5) throw new Error(`ATOMS: рядок має 5 полів: ${r}`);
        const [id_s, Z_s, x_s, y_s, z_s] = parts;
        const idx = parseInt(id_s, 10);
        const Z = parseInt(Z_s, 10);
        const x = parseFloat(x_s), y = parseFloat(y_s), z = parseFloat(z_s);
        if ([idx, Z, x, y, z].some(v => Number.isNaN(v))) throw new Error(`ATOMS: некоректні типи: ${r}`);
        if (idx < 0) throw new Error(`ATOMS: id має бути ≥0 (отримано ${idx}).`);
        if (seen.has(idx)) throw new Error(`ATOMS: дубльований id=${idx}.`);
        seen.add(idx);
        atoms.push([idx, { Z, x, y, z }]);
    }
    if (!atoms.length) throw new Error("ATOMS: порожній список атомів.");
    atoms.sort((a, b) => a[0] - b[0]);
    atoms.forEach((pair, expect) => {
        if (pair[0] !== expect) throw new Error(`ATOMS: очікувався id=${expect}, отримано id=${pair[0]}.`);
    });
    const atomsOut = atoms.map(p => p[1]);
    const N = atomsOut.length;

    const bonds = [];
    if (i_bonds >= 0) {
        if (i_bonds + 1 >= lines.length) throw new Error("Очікується заголовок 'i,j,order' після BONDS.");
        if (lines[i_bonds + 1] !== "i,j,order") throw new Error("Невірний заголовок BONDS.");

        const bond_rows = lines.slice(i_bonds + 2).filter(Boolean);
        for (const r of bond_rows) {
            const parts = r.split(",");
            if (parts.length !== 3) throw new Error(`BONDS: рядок має 3 поля: ${r}`);
            const i = parseInt(parts[0], 10), j = parseInt(parts[1], 10), o = parseInt(parts[2], 10);
            if ([i, j, o].some(Number.isNaN)) throw new Error(`BONDS: некоректні типи: ${r}`);
            if (!(i >= 0 && i < N && j >= 0 && j < N)) throw new Error(`BONDS: індекси поза межами [0..${N - 1}]: ${r}`);
            if (i === j) throw new Error(`BONDS: i та j не можуть збігатися: ${r}`);
            if (![1, 2, 3, 4].includes(o)) throw new Error(`BONDS: order має бути 1|2|3|4: ${r}`);
            bonds.push([i, j, o]);
        }
    }
    return [atomsOut, bonds];
}
