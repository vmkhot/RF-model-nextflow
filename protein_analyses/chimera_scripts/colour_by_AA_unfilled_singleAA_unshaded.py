from chimerax.core.commands import run
import csv, sys

AA_THREE_LETTER = {
    "R": "ARG", "H": "HIS", "K": "LYS",
    "D": "ASP", "E": "GLU",
    "S": "SER", "T": "THR", "N": "ASN", "Q": "GLN",
    "A": "ALA", "V": "VAL", "I": "ILE", "L": "LEU", "M": "MET",
    "F": "PHE", "Y": "TYR", "W": "TRP",
    "C": "CYS", "U": "SEC", "G": "GLY", "P": "PRO",
}

AA_GROUPS = {
    "R": "#ef4806", "H": "#ef4806", "K": "#ef4806",
    "D": "#234abb", "E": "#234abb",
    "S": "#1eae83", "T": "#1eae83", "N": "#1eae83", "Q": "#1eae83",
    "A": "#f4b713", "V": "#f4b713", "I": "#f4b713", "L": "#f4b713", "M": "#f4b713",
    "F": "#bf21a2", "Y": "#bf21a2", "W": "#bf21a2",
    "C": "#4a4444", "U": "#4a4444", "G": "#4a4444", "P": "#4a4444",
}

# -------------------------
# Main logic
# -------------------------
def run_script(session, tsv_file=None, model="#1", biome=None, aa_filter=None, color_override=None):

    # -------------------------
    # Native mode
    # -------------------------
    if biome is None or biome == "native":
        print(">>> Coloring residues (native mode)")

        for aa, group_color in AA_GROUPS.items():
            if aa_filter and aa != aa_filter:
                continue

            three = AA_THREE_LETTER.get(aa)
            if three is None:
                continue

            hex_color = color_override if color_override else group_color

            cmd = f"color {model}:{three} {hex_color}"
            print(f">>> Running: {cmd}")
            try:
                run(session, cmd)
            except Exception as e:
                print(f"WARNING: failed for {three}: {e}")

        print(">>> Coloring complete!")
        return

    # -------------------------
    # TSV mode
    # -------------------------
    if tsv_file is None:
        print("ERROR: Must provide tsv_file when using TSV mode")
        return

    print(f">>> Loading: {tsv_file}, filtering biome='{biome}'")

    pos_to_data = {}

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

                if aa_filter and aa != aa_filter:
                    continue

                if aa and aa in AA_GROUPS:
                    pos_to_data[p] = aa

            except Exception:
                pass

    if not pos_to_data:
        print(f"ERROR: No data loaded for biome='{biome}'")
        return

    print(f">>> Loaded {len(pos_to_data)} positions")

    # -------------------------
    # Apply coloring (flat, no importance-based shading)
    # -------------------------
    for pos, aa in sorted(pos_to_data.items()):
        hex_color = color_override if color_override else AA_GROUPS[aa]

        cmd = f"color {model}:{pos} {hex_color}"
        print(f">>> Running: {cmd} ({aa})")
        run(session, cmd)

    print(">>> Coloring complete!")

# -------------------------
# CLI handling
# -------------------------
mode_or_file  = sys.argv[1]
model_id      = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome         = sys.argv[3] if len(sys.argv) > 3 else None
aa_filter     = sys.argv[4] if len(sys.argv) > 4 else None
color_override = sys.argv[5] if len(sys.argv) > 5 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native", aa_filter=aa_filter, color_override=color_override)
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome, aa_filter=aa_filter, color_override=color_override)



# native mode, default AA-group colors
# python script.py native #1

# native mode, only color Phenylalanine residues, custom color
# python script.py native #1 native F "#ff00ff"

# TSV mode, ocean biome, all AAs, default group colors (no shading)
# python script.py mydata.tsv #1 ocean

# TSV mode, ocean biome, only Arg residues, custom color
# python script.py mydata.tsv #1 ocean R "#ff00ff"