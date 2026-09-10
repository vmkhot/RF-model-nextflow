from chimerax.core.commands import run
from chimerax.atomic import all_atomic_structures
import csv, sys

AA_THREE_TO_ONE = {
    "ALA": "A", "VAL": "V", "ILE": "I", "LEU": "L", "MET": "M",
    "PHE": "F", "TYR": "Y", "TRP": "W", "SER": "S", "THR": "T",
    "ASN": "N", "GLN": "Q", "CYS": "C", "SEC": "U", "GLY": "G",
    "PRO": "P", "ARG": "R", "HIS": "H", "LYS": "K", "ASP": "D",
    "GLU": "E",
}

AA_THREE_LETTER = {v: k for k, v in AA_THREE_TO_ONE.items()}

AA_GROUPS = {
    "R": "#ef4806", "H": "#ef4806", "K": "#ef4806",
    "D": "#234abb", "E": "#234abb",
    "S": "#1eae83", "T": "#1eae83", "N": "#1eae83", "Q": "#1eae83",
    "A": "#f4b713", "V": "#f4b713", "I": "#f4b713", "L": "#f4b713", "M": "#f4b713",
    "F": "#bf21a2", "Y": "#bf21a2", "W": "#bf21a2",
    "C": "#4a4444", "U": "#4a4444", "G": "#4a4444", "P": "#4a4444",
}

def get_structure_residues(session, model):
    """Returns dict of {seq_pos: one_letter_aa} from the loaded structure."""
    structures = all_atomic_structures(session)
    if not structures:
        print("ERROR: No structures loaded in session")
        return {}

    # Match model id if specified
    model_num = model.replace("#", "")
    target = None
    for s in structures:
        if str(s.id[0]) == model_num:
            target = s
            break
    if target is None:
        target = structures[0]
        print(f"WARNING: Model {model} not found, using first loaded structure")

    residue_map = {}
    for res in target.residues:
        one = AA_THREE_TO_ONE.get(res.name.upper())
        if one:
            residue_map[res.number] = one

    print(f">>> Found {len(residue_map)} residues in structure")
    return residue_map


def run_script(session, tsv_file=None, model="#1", biome=None):
    if biome is None or biome == "native":
        print(">>> Coloring all residues by amino acid group (native mode)")
        for aa, hex_color in AA_GROUPS.items():
            three = AA_THREE_LETTER.get(aa)
            if three is None:
                continue
            cmd = f"color {model}:{three} {hex_color}"
            print(f">>> Running: {cmd}")
            try:
                run(session, cmd)
            except Exception as e:
                print(f"    WARNING: failed for {three}: {e}")
        print(">>> Coloring complete!")
        return

    # TSV mode
    if tsv_file is None:
        print("ERROR: Must provide tsv_file=/path/to/file when using lake/ocean mode")
        return

    print(f">>> Loading: {tsv_file}, filtering by biome='{biome}'")

    # 1. Get all residues from the loaded structure
    structure_residues = get_structure_residues(session, model)
    if not structure_residues:
        return

    # 2. Read TSV positions filtered by biome
    tsv_pos_to_aa = {}
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                parts = row["feature"].split("_")
                p = int(float(parts[0]) + 1)
                aa = parts[1].strip().upper() if len(parts) > 1 else None
                row_biome = row["biome_predictor"].strip().lower()

                if biome.lower() != "native" and row_biome != biome.lower():
                    continue

                if aa and aa in AA_GROUPS:
                    tsv_pos_to_aa[p] = aa
            except:
                pass

    print(f">>> {len(tsv_pos_to_aa)} positions found in TSV for biome='{biome}'")

    # 3. Build final coloring map:
    #    - TSV positions → use AA from TSV feature
    #    - All other structure positions → use native AA from structure
    final_pos_to_aa = {}
    for pos, native_aa in structure_residues.items():
        if pos in tsv_pos_to_aa:
            final_pos_to_aa[pos] = tsv_pos_to_aa[pos]  # TSV takes priority
        else:
            final_pos_to_aa[pos] = native_aa            # fall back to native

    print(f">>> Coloring {len(final_pos_to_aa)} total positions")
    print(f">>>   {len(tsv_pos_to_aa)} from TSV, "
          f"{len(final_pos_to_aa) - len(tsv_pos_to_aa)} from native structure")

    # 4. Color each position
    for pos, aa in sorted(final_pos_to_aa.items()):
        hex_color = AA_GROUPS.get(aa)
        if hex_color is None:
            continue
        cmd = f"color {model}:{pos} {hex_color}"
        run(session, cmd)

    print(">>> Coloring complete!")


mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome)