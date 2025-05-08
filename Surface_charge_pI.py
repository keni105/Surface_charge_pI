#!/usr/bin/env python3
"""
Rank protein models by net surface-charge and compute pI from structure (PROPKA)
or sequence if PROPKA is unavailable.  Supports both .cif and .pdb inputs,
uses APBS when available, and falls back to a simple residue-based charge.
"""
import os
import subprocess
import glob
from collections import Counter

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# three-letter → one-letter map for standard amino acids
_three_to_one = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

def convert_cif_to_pdb(cif_file, output_dir="converted_pdb"):
    """Convert an mmCIF file to PDB format using Biopython."""
    try:
        os.makedirs(output_dir, exist_ok=True)
        parser = MMCIFParser(QUIET=True)
        struct_id = os.path.splitext(os.path.basename(cif_file))[0]
        structure = parser.get_structure(struct_id, cif_file)
        pdb_path = os.path.join(output_dir, f"{struct_id}.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        return pdb_path
    except Exception as e:
        print(f"✗  CIF→PDB failed for {cif_file}: {e}")
        return None

def fallback_charge(pdb_file):
    """Residue-based net charge at pH≈7."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)
    counts = Counter(res.resname for res in structure.get_residues() if res.id[0]==' ')
    return (
        counts.get('ARG',0)* 1.0 +
        counts.get('LYS',0)* 1.0 +
        counts.get('ASP',0)*-1.0 +
        counts.get('GLU',0)*-1.0 +
        counts.get('HIS',0)* 0.1
    )

def calculate_surface_charge(pdb_file):
    """Try APBS to get net charge; else residue-based fallback."""
    base = os.path.splitext(os.path.basename(pdb_file))[0]
    pqr = f"{base}.pqr"
    # 1) PDB→PQR
    try:
        subprocess.run(
            ["pdb2pqr", "--ff=PARSE", pdb_file, pqr],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except Exception:
        return fallback_charge(pdb_file), 'residue-count'
    # 2) Minimal APBS input
    in_file = f"{base}.in"
    with open(in_file, 'w') as inp:
        inp.write(f"""
read
    mol pqr {pqr}
end
energy
    calculate total charge
end
quit
""")
    # 3) Run APBS
    try:
        out = subprocess.run(
            ["apbs", in_file],
            capture_output=True, text=True, check=True
        ).stdout
        for line in out.splitlines():
            if "Net charge" in line:
                charge = float(line.split()[-2])
                return charge, 'APBS'
        raise RuntimeError
    except Exception:
        print(f"→ Using fallback residue-based estimate for {pdb_file}")
        return fallback_charge(pdb_file), 'residue-count'

def calculate_pI(pdb_file):
    """
    First try PROPKA (structure‐based).  If that fails, use sequence‐based ProtParam.
    """
    base = os.path.splitext(os.path.basename(pdb_file))[0]
    # try PROPKA
    try:
        # run propka3 → writes base.pka
        subprocess.run(["propka3", pdb_file], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pka_file = f"{base}.pka"
        if os.path.isfile(pka_file):
            with open(pka_file) as f:
                for L in f:
                    if "Isoelectric point" in L or "pI" in L:
                        # line like "    Isoelectric point (pI)   5.23"
                        parts = L.strip().split()
                        val = [p for p in parts if p.replace('.','',1).isdigit()]
                        if val:
                            return float(val[-1]), 'PROPKA'
    except Exception:
        pass

    # fallback: sequence-based
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)
    seq = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0]==' ':
                    aa = _three_to_one.get(res.resname)
                    if aa:
                        seq.append(aa)
    seq = "".join(seq)
    if not seq:
        return float('nan'), 'sequence'
    pi = ProteinAnalysis(seq).isoelectric_point()
    return pi, 'sequence'

def rank_proteins(files):
    """
    Process .cif/.pdb list → [(name, charge, charge_method, pI, pI_method)]
    sorted by |charge| desc.
    """
    out=[]
    for fn in files:
        print(f"→ Processing {fn}")
        ext = os.path.splitext(fn)[1].lower()
        if ext=='.cif':
            pdb = convert_cif_to_pdb(fn)
            if not pdb: continue
        elif ext=='.pdb':
            pdb=fn
        else:
            print(f"✗  Skipping {fn}")
            continue
        q, cm = calculate_surface_charge(pdb)
        pi, pim = calculate_pI(pdb)
        out.append((os.path.basename(fn), q, cm, pi, pim))
    return sorted(out, key=lambda x: abs(x[1]), reverse=True)

if __name__=='__main__':
    files = sorted(glob.glob("*.cif") + glob.glob("*.pdb"))
    if not files:
        print("No .cif/.pdb found."); exit(1)
    ranked = rank_proteins(files)
    print("\nRanked Proteins by surface-charge + pI:")
    for i,(name,q,cm,pi,pim) in enumerate(ranked,1):
        print(f"{i}. {name:<40} Q={q:7.2f} ({cm})  pI={pi:5.2f} ({pim})")
