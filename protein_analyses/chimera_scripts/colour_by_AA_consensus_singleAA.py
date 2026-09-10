from chimerax.core.commands import run
from chimerax.atomic import all_atomic_structures
from Bio import SeqIO
import sys

AA_GROUPS = {
    "R": "#ef4806", "H": "#ef4806", "K": "#ef4806",
    "D": "#234abb", "E": "#234abb",
    "S": "#1eae83", "T": "#1eae83", "N": "#1eae83", "Q": "#1eae83",
    "A": "#f4b713", "V": "#f4b713", "I": "#f4b713", "L": "#f4b713", "M": "#f4b713",
    "F": "#bf21a2", "Y": "#bf21a2", "W": "#bf21a2",
    "C": "#4a4444", "U": "#4a4444", "G": "#4a4444", "P": "#4a4444",
    "X": "#e5e5e5",
}

def get_structure_residues(session, model):
    structures = all_atomic_structures(session)
    if not structures:
        print("ERROR: No structures loaded in session")
        return {}

    model_num = model.replace("#", "")
    target = None
    for s in structures:
        if str(s.id[0]) == model_num:
            target = s
            break
    if target is None:
        target = structures[0]
        print(f"WARNING: Model {model} not found, using first loaded structure")

    residue_list = list(target.residues)
    print(f">>> Found {len(residue_list)} residues in structure")
    return residue_list


def parse_consensus(fasta_file, biome):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if biome.lower() in rec.id.lower():
            print(f">>> Matched record: {rec.id}")
            print(f">>> Consensus length: {len(rec.seq)}")
            return str(rec.seq).upper()
    return None


def run_script(session, fasta_file=None, model="#1", biome="Lake", target_aa=None):
    if fasta_file is None:
        print("ERROR: Must provide fasta_file=/path/to/consensus.fasta")
        return

    if target_aa:
        target_aa = target_aa.upper()
        if len(target_aa) != 1:
            print("ERROR: target_aa must be a single letter")
            return
        print(f">>> Filtering to amino acid: {target_aa}")

    print(f">>> Loading consensus from: {fasta_file}")
    print(f">>> Biome: {biome}")

    consensus = parse_consensus(fasta_file, biome)
    if consensus is None:
        print(f"ERROR: No consensus sequence found for biome '{biome}'")
        return

    residue_list = get_structure_residues(session, model)
    if not residue_list:
        return

    if len(consensus) != len(residue_list):
        print(f"WARNING: Consensus length ({len(consensus)}) != "
              f"structure residues ({len(residue_list)})")
        print(f"         Will color up to min({len(consensus)}, {len(residue_list)}) positions")

    n = min(len(consensus), len(residue_list))

    print(f">>> Coloring {n} residues...")
    unknown_aas = set()

    for i in range(n):
        aa = consensus[i]

        if target_aa is not None and aa != target_aa:
            continue

        residue = residue_list[i]
        res_num = residue.number
        hex_color = AA_GROUPS.get(aa)

        if hex_color is None:
            unknown_aas.add(aa)
            hex_color = "#e5e5e5"

        cmd = f"color {model}:{res_num} {hex_color}"
        run(session, cmd)

    if unknown_aas:
        print(f"WARNING: Unknown AA characters treated as grey: {unknown_aas}")

    print("\n>>> Color legend:")
    groups = {
        "Positive charge (R,H,K)" : "#ef4806",
        "Negative charge (D,E)"   : "#234abb",
        "Polar uncharged (S,T,N,Q)": "#1eae83",
        "Hydrophobic (A,V,I,L,M)" : "#f4b713",
        "Aromatic (F,Y,W)"        : "#bf21a2",
        "Other (C,U,G,P)"         : "#4a4444",
        "Ambiguous (X)"           : "#e5e5e5",
    }
    for label, color in groups.items():
        print(f"    {color}  {label}")

    print("\n>>> Coloring complete!")


# --- Entry point ---
# Usage:
# chimerax --script "script.py consensus.fasta #1 Lake"
# chimerax --script "script.py consensus.fasta #1 Lake A"

fasta_file = sys.argv[1]
model_id   = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome      = sys.argv[3] if len(sys.argv) > 3 else "Lake"
target_aa  = sys.argv[4] if len(sys.argv) > 4 else None

run_script(session, fasta_file=fasta_file, model=model_id, biome=biome, target_aa=target_aa)