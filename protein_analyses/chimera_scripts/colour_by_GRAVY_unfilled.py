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

# Kyte & Doolittle (1982) hydropathy values
KD_HYDROPATHY = {
    "I": 4.5,  "V": 4.2,  "L": 3.8,  "F": 2.8,  "C": 2.5,
    "M": 1.9,  "A": 1.8,  "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
    "U": 2.5,  # selenocysteine, treated like cysteine
}

GRAVY_MIN = -4.5
GRAVY_MAX = 4.5

# ColorBrewer BrBG-11 diverging palette (low=brown/hydrophilic,
# mid=white, high=teal-blue/hydrophobic), as (R, G, B) stops evenly
# spaced across the GRAVY range.
BRBG_STOPS = [
    (0x00, 0x3c, 0x30),
    (0x01, 0x66, 0x5e),
    (0x35, 0x97, 0x8f),
    (0x80, 0xcd, 0xc1),
    (0xc7, 0xea, 0xe5),
    (0xf5, 0xf5, 0xf5),
    (0xf6, 0xe8, 0xc3),
    (0xdf, 0xc2, 0x7d),
    (0xbf, 0x81, 0x2d),
    (0x8c, 0x51, 0x0a),
    (0x54, 0x30, 0x05),
]

def gravy_to_hex(value):
    t = (value - GRAVY_MIN) / (GRAVY_MAX - GRAVY_MIN)
    t = max(0.0, min(1.0, t))
    n = len(BRBG_STOPS) - 1
    scaled = t * n
    i = min(int(scaled), n - 1)
    frac = scaled - i
    c0, c1 = BRBG_STOPS[i], BRBG_STOPS[i + 1]
    r = round(c0[0] + frac * (c1[0] - c0[0]))
    g = round(c0[1] + frac * (c1[1] - c0[1]))
    b = round(c0[2] + frac * (c1[2] - c0[2]))
    return f"#{r:02x}{g:02x}{b:02x}"

def aa_color(aa):
    value = KD_HYDROPATHY.get(aa)
    if value is None:
        return None
    return gravy_to_hex(value)

# Name of the column in the TSV holding the importance/gini score used to
# rank rows so we know which one is "top" for a repeated position.
GINI_COLUMN = "feature_importance_vals"

def run_script(session, tsv_file=None, model="#1", biome=None):
    if biome is None or biome == "native":
        print(">>> Coloring all residues by GRAVY (Kyte-Doolittle) gradient (native mode)")
        for aa, three in AA_THREE_LETTER.items():
            hex_color = aa_color(aa)
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

    rows = []
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
                if not aa or aa_color(aa) is None:
                    continue
                gini = float(row[GINI_COLUMN])
                rows.append((gini, p, aa))
            except Exception:
                pass

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
        hex_color = aa_color(aa)
        cmd = f"color {model}:{pos} {hex_color}"
        print(f">>> Running: {cmd} ({aa}, GRAVY={KD_HYDROPATHY[aa]})")
        run(session, cmd)

    print(">>> Coloring complete!")

mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome)