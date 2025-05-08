# SurfaceCharge_pI

A Python utility to **rank protein models** by their net surface-charge and compute their isoelectric point (pI) from structure (via PROPKA) or sequence fallback. It accepts **both** mmCIF (`.cif`) and PDB (`.pdb`) files.

---

## 🔧 Features

1. **Surface-charge calculation**
   - Attempts a full electrostatics solve with **APBS** (via `pdb2pqr → APBS`).
   - Falls back to a **simple residue-based estimate** if APBS or `pdb2pqr` fails.

2. **pI estimation**
   - **PROPKA** structure-based pI (via `propka3 <pdb>` → `.pka` file).
   - Falls back to **sequence-based** pI using Biopython's `ProteinAnalysis`.

3. **Input flexibility**
   - Automatically converts `.cif` → `.pdb` (Biopython).
   - Supports native `.pdb` inputs unchanged.

4. **Automated ranking & reporting**
   - Sorts models by **absolute surface-charge** descending.
   - Prints a table: filename, charge (method), pI (method).

---

## 📦 Requirements

- **Python 3.7+**
- **Biopython**
- **pdb2pqr** (PARSE force field)
- **APBS** (≥3.0)
- **propka3** (PROPKA 3.x)

You can install core tools via `conda` or `brew`:

```bash
# Create env
conda create -n surfpi python=3.9 -y
conda activate surfpi

# Python deps
pip install biopython

# pdb2pqr + APBS
conda install -c conda-forge pdb2pqr apbs -y
# or on macOS:
brew install apbs pdb2pqr

# PROPKA
# If available on conda-forge:
# conda install -c conda-forge propka3
# Otherwise, clone & install manually:
# git clone https://github.com/jensengroup/propka.git
# cd propka/PROPKA-3.1 && python setup.py install
```

Ensure all binaries are on your `PATH`:
```bash
which pdb2pqr apbs propka3
```

---

## 🚀 Usage

1. **Place** your `.cif` and/or `.pdb` files in a directory.
2. **Run**:
   ```bash
   python Surface_charge_pI.py
   ```
3. **Result** printed to console, e.g.:
   ```text
   Ranked Proteins by surface-charge + pI:
   1. modelA.cif Q=-17.00 (APBS)  pI=5.10 (PROPKA)
   2. modelB.pdb Q=-14.00 (residue-count)  pI=6.20 (sequence)
   ...
   ```

Each line shows:
- **Q** = net surface-charge (in e)
- **(APBS)** or **(residue-count)** method used
- **pI** and whether **(PROPKA)** or **(sequence)** was used

---

## ⚠️ Drawbacks & Caveats

- **APBS grid parameters** are minimal and may not suit very large or complex structures.
- **PROPKA** structural pI depends on having correctly folded/complete PDBs.
- **Fallback estimates** ignore local pKa shifts and detailed solvation.
- **Performance**: APBS and PROPKA runs can be slow for multiple, large proteins.
- **External dependencies**: Binaries (`pdb2pqr`, `apbs`, `propka3`) must be installed and on `PATH`.

---

## 📄 License

MIT © Keni Vidilaseris
