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


def hex_to_rgb(hex_color):
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))


def blend_with_white(hex_color, t):
    """t=0 → white, t=1 → full color."""
    r, g, b = hex_to_rgb(hex_color)
    return rgb_to_hex(1.0 + t*(r-1.0), 1.0 + t*(g-1.0), 1.0 + t*(b-1.0))


def normalize(vals_dict):
    vals  = list(vals_dict.values())
    vmin, vmax = min(vals), max(vals)
    denom = vmax - vmin + 1e-12
    return {k: (v - vmin) / denom for k, v in vals_dict.items()}


def run_script(session, tsv_file=None, model="#1", biome=None,
               min_shade=0.15, max_shade=1.0):

    # ── NATIVE MODE ──────────────────────────────────────────────────────────
    if biome is None or biome == "native":
        print(">>> Native mode: flat AA group colors")
        for aa, hex_color in AA_GROUPS.items():
            three = AA_THREE_LETTER.get(aa)
            if not three:
                continue
            try:
                run(session, f"color {model}:{three} {hex_color}")
            except Exception as e:
                print(f"    WARNING: {three}: {e}")
        print(">>> Coloring complete!")
        return

    # ── TSV MODE ─────────────────────────────────────────────────────────────
    if tsv_file is None:
        print("ERROR: Must provide tsv_file when using lake/ocean mode")
        return

    print(f">>> Loading: {tsv_file}  biome='{biome}'")

    pos_to_aa  = {}
    pos_to_imp = {}

    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                parts     = row["feature"].split("_")
                p         = int(float(parts[0]) + 1)
                aa        = parts[1].strip().upper() if len(parts) > 1 else None
                row_biome = row["biome_predictor"].strip().lower()
                imp       = float(row["feature_importance_vals"])

                if biome.lower() != "native" and row_biome != biome.lower():
                    continue
                if not aa or aa not in AA_GROUPS:
                    continue

                if p not in pos_to_imp or imp > pos_to_imp[p]:
                    pos_to_aa[p]  = aa
                    pos_to_imp[p] = imp
            except:
                pass

    if not pos_to_aa:
        print(f"ERROR: No data loaded for biome='{biome}'")
        return

    print(f">>> Loaded {len(pos_to_aa)} positions — coloring now...")

    norm_imp = normalize(pos_to_imp)

    # Color only TSV positions — everything else stays ChimeraX default
    for pos, aa in sorted(pos_to_aa.items()):
        t            = min_shade + norm_imp[pos] * (max_shade - min_shade)
        shaded_color = blend_with_white(AA_GROUPS[aa], t)
        run(session, f"color {model}:{pos} {shaded_color}")

    print(f">>> Coloring complete! ({len(pos_to_aa)} residues colored)")


# ── ENTRY POINT ───────────────────────────────────────────────────────────────
# Usage:
#   Native (flat AA colors, no importance):
#     chimerax --script "script.py native #1"
#
#   TSV + biome (AA colors shaded by importance):
#     chimerax --script "script.py data.tsv #1 lake"
#     chimerax --script "script.py data.tsv #1 ocean"
#
#   Optional: control shading range (default 0.15 1.0)
#     chimerax --script "script.py data.tsv #1 lake 0.2 1.0"

# for shading: 0.15 = minimal saturation (very light, low importance)
# for shading: 1 = full saturation (full, high importance)

mode_or_file = sys.argv[1]
model_id     = sys.argv[2] if len(sys.argv) > 2 else "#1"
biome        = sys.argv[3] if len(sys.argv) > 3 else None
min_shade    = float(sys.argv[4]) if len(sys.argv) > 4 else 0.15
max_shade    = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0

if mode_or_file == "native":
    run_script(session, model=model_id, biome="native")
else:
    run_script(session, tsv_file=mode_or_file, model=model_id,
               biome=biome, min_shade=min_shade, max_shade=max_shade)