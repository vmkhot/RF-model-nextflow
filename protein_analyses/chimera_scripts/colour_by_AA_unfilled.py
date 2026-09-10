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

# Name of the column in the TSV holding the importance/gini score used to
# rank rows so we know which one is "top" for a repeated position.
GINI_COLUMN = "feature_importance_vals"

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

    rows = []
    n_total = 0
    n_biome_mismatch = 0
    n_bad_aa = 0
    n_errors = 0
    error_examples = []
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        print(f">>> Detected columns: {reader.fieldnames}")
        for row in reader:
            n_total += 1
            try:
                parts = row["feature"].split("_")
                p = int(float(parts[0]) + 1)
                aa = parts[1].strip().upper() if len(parts) > 1 else None
                row_biome = row["biome_predictor"].strip().lower()
                if biome.lower() != "native" and row_biome != biome.lower():
                    n_biome_mismatch += 1
                    continue
                if not aa or aa not in AA_GROUPS:
                    n_bad_aa += 1
                    continue
                gini = float(row[GINI_COLUMN])
                rows.append((gini, p, aa))
            except Exception as e:
                n_errors += 1
                if len(error_examples) < 5:
                    error_examples.append((row, str(e)))

    print(f">>> Rows read: {n_total}, kept: {len(rows)}, "
          f"biome mismatch: {n_biome_mismatch}, bad/missing aa: {n_bad_aa}, "
          f"errors: {n_errors}")
    if error_examples:
        print(">>> Example errors (up to 5):")
        for row, err in error_examples:
            print(f"    row={row} -> {err}")

    if not rows:
        print(f"ERROR: No data loaded for biome='{biome}'")
        return

    # Sort by gini descending so the highest-importance row for each
    # position is encountered first.
    rows.sort(key=lambda r: r[0], reverse=True)

    pos_to_aa = {}
    for gini, p, aa in rows:
        if p in pos_to_aa:
            # Already assigned from a higher (or equal, earlier) gini row -
            # skip this lower-ranked duplicate for the same position.
            continue
        pos_to_aa[p] = aa

    print(f">>> Loaded {len(pos_to_aa)} positions (top gini per position)")

    for pos, aa in sorted(pos_to_aa.items()):
        hex_color = AA_GROUPS[aa]
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