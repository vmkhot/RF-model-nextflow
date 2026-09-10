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

# Hydrophobic (nonpolar) residues -> yellow
# Hydrophilic (polar/charged) residues -> blue
HYDROPHOBIC = {"A", "V", "L", "I", "P", "F", "M", "W", "G", "C", "U"}
HYDROPHILIC = {"R", "H", "K", "D", "E", "S", "T", "N", "Q", "Y"}

HYDROPHOBIC_COLOR = "#f4d613"  # yellow
HYDROPHILIC_COLOR = "#234abb"  # blue

def classify_color(aa):
    if aa in HYDROPHOBIC:
        return HYDROPHOBIC_COLOR
    elif aa in HYDROPHILIC:
        return HYDROPHILIC_COLOR
    return None

def run_script(session, tsv_file=None, model="#1", biome=None):
    if biome is None or biome == "native":
        print(">>> Coloring all residues by hydrophobicity (native mode)")
        for aa, three in AA_THREE_LETTER.items():
            hex_color = classify_color(aa)
            if hex_color is None:
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
    pos_to_aa = {}
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
                if aa and (aa in HYDROPHOBIC or aa in HYDROPHILIC):
                    pos_to_aa[p] = aa
            except:
                pass

    if not pos_to_aa:
        print(f"ERROR: No data loaded for biome='{biome}'")
        return

    print(f">>> Loaded {len(pos_to_aa)} positions")
    for pos, aa in sorted(pos_to_aa.items()):
        hex_color = classify_color(aa)
        cmd = f"color {model}:{pos} {hex_color}"
        print(f">>> Running: {cmd} ({aa})")
        run(session, cmd)
    print(">>> Coloring complete!")

mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome)