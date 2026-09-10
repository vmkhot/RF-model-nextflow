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

# Approximate side-chain charge at physiological pH (~7.4).
# Acidic residues (D, E) are negative; basic residues (K, R) are positive;
# histidine is mostly neutral but carries a small partial positive charge
# (pKa ~6, so a minority is protonated at pH 7.4); everything else is
# treated as uncharged for this coloring scheme.
CHARGE = {
    "D": -1.0, "E": -1.0,          # acidic
    "K": 1.0,  "R": 1.0,           # basic
    "H": 0.1,                      # weakly/partially basic
    "S": 0.0, "T": 0.0, "N": 0.0, "Q": 0.0,
    "A": 0.0, "V": 0.0, "I": 0.0, "L": 0.0, "M": 0.0,
    "F": 0.0, "Y": 0.0, "W": 0.0,
    "C": 0.0, "U": 0.0, "G": 0.0, "P": 0.0,
}

CHARGE_MIN = -1.0
CHARGE_MAX = 1.0

# ColorBrewer RdBu-11 diverging palette (low=dark red/negative,
# mid=white/neutral, high=dark blue/positive), as (R, G, B) stops evenly
# spaced across the charge range. This matches the common electrostatics
# convention (e.g. APBS/PyMOL): red = negative, blue = positive.
RDBU_STOPS = [
    (0x67, 0x00, 0x1f),
    (0xb2, 0x18, 0x2b),
    (0xd6, 0x60, 0x4d),
    (0xf4, 0xa5, 0x82),
    (0xfd, 0xdb, 0xc7),
    (0xf7, 0xf7, 0xf7),
    (0xd1, 0xe5, 0xf0),
    (0x92, 0xc5, 0xde),
    (0x43, 0x93, 0xc3),
    (0x21, 0x66, 0xac),
    (0x05, 0x30, 0x61),
]

def charge_to_hex(value):
    t = (value - CHARGE_MIN) / (CHARGE_MAX - CHARGE_MIN)
    t = max(0.0, min(1.0, t))
    n = len(RDBU_STOPS) - 1
    scaled = t * n
    i = min(int(scaled), n - 1)
    frac = scaled - i
    c0, c1 = RDBU_STOPS[i], RDBU_STOPS[i + 1]
    r = round(c0[0] + frac * (c1[0] - c0[0]))
    g = round(c0[1] + frac * (c1[1] - c0[1]))
    b = round(c0[2] + frac * (c1[2] - c0[2]))
    return f"#{r:02x}{g:02x}{b:02x}"

def aa_color(aa):
    value = CHARGE.get(aa)
    if value is None:
        return None
    return charge_to_hex(value)

# Name of the column in the TSV holding the importance/gini score used to
# rank rows so we know which one is "top" for a repeated position.
GINI_COLUMN = "feature_importance_vals"

def run_script(session, tsv_file=None, model="#1", biome=None):
    if biome is None or biome == "native":
        print(">>> Coloring all residues by charge (electrostatics) gradient (native mode)")
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
        print(f">>> Running: {cmd} ({aa}, charge={CHARGE[aa]})")
        run(session, cmd)

    print(">>> Coloring complete!")

mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id, biome=biome)