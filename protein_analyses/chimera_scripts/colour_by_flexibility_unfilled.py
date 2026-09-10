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

# Vihinen, Torkkila & Riikonen (1994) normalized average flexibility
# parameters (B-value derived). Proteins 19(2):141-9.
VIHINEN_FLEX = {
    "A": 0.984, "C": 0.906, "E": 1.094, "D": 1.068,
    "G": 1.031, "F": 0.915, "I": 0.927, "H": 0.950,
    "K": 1.102, "M": 0.952, "L": 0.935, "N": 1.048,
    "Q": 1.037, "P": 1.049, "S": 1.046, "R": 1.008,
    "T": 0.997, "W": 0.904, "V": 0.931, "Y": 0.929,
    "U": 0.906,  # selenocysteine, treated like cysteine
}

# Bounds of the Vihinen scale (Trp = min, Lys = max), used to normalize
# values onto the [0, 1] range for color interpolation.
FLEX_MIN = 0.904
FLEX_MAX = 1.102

# ColorBrewer PiYG-11 diverging palette (low=pink/less flexible,
# mid=white, high=green/more flexible), as (R, G, B) stops evenly
# spaced across the flexibility range.
PIYG_STOPS = [
    (0x8e, 0x01, 0x52),
    (0xc5, 0x1b, 0x7d),
    (0xde, 0x77, 0xae),
    (0xf1, 0xb6, 0xda),
    (0xfd, 0xe0, 0xef),
    (0xf7, 0xf7, 0xf7),
    (0xe6, 0xf5, 0xd0),
    (0xb8, 0xe1, 0x86),
    (0x7f, 0xbc, 0x41),
    (0x4d, 0x92, 0x21),
    (0x27, 0x64, 0x19),
]

def flex_to_hex(value):
    t = (value - FLEX_MIN) / (FLEX_MAX - FLEX_MIN)
    t = max(0.0, min(1.0, t))
    n = len(PIYG_STOPS) - 1
    scaled = t * n
    i = min(int(scaled), n - 1)
    frac = scaled - i
    c0, c1 = PIYG_STOPS[i], PIYG_STOPS[i + 1]
    r = round(c0[0] + frac * (c1[0] - c0[0]))
    g = round(c0[1] + frac * (c1[1] - c0[1]))
    b = round(c0[2] + frac * (c1[2] - c0[2]))
    return f"#{r:02x}{g:02x}{b:02x}"

def aa_color(aa):
    value = VIHINEN_FLEX.get(aa)
    if value is None:
        return None
    return flex_to_hex(value)

# Name of the column in the TSV holding the importance/gini score used to
# rank rows so we know which one is "top" for a repeated position.
GINI_COLUMN = "feature_importance_vals"

def run_script(session, tsv_file=None, model="#1", biome=None):
    if biome is None or biome == "native":
        print(">>> Coloring all residues by Vihinen flexibility gradient (native mode)")
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
        print(f">>> Running: {cmd} ({aa}, Flex={VIHINEN_FLEX[aa]})")
        run(session, cmd)

    print(">>> Coloring complete!")

mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome)