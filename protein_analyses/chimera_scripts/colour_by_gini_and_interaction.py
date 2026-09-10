from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import csv, sys


def get_all_residue_numbers(session, model_spec):
    """
    Return a sorted set of all residue numbers (integers) that exist
    in the model matching model_spec, across all chains.
    """
    all_positions = set()
    for m in session.models:
        if isinstance(m, AtomicStructure):
            if model_spec.lstrip("#") == "" or str(m.id_string) == model_spec.lstrip("#"):
                for r in m.residues:
                    all_positions.add(r.number)
    return all_positions


def run_script(session, tsv_file=None, model="#1",
               interaction_type="all", palette_color="a70088"):

    if tsv_file is None:
        print("ERROR: Must provide tsv_file=/path/to/file")
        return

    print(">>> Loading:", tsv_file)

    # 1. Read every row, tracking:
    #    - full_vals: ALL feature_importance_vals_A across the whole file,
    #      keyed by position (used only for normalization range)
    #    - pos_to_val: feature_importance_vals_A for positions that match
    #      the requested interaction_type (used for what gets colored)
    full_vals = {}
    pos_to_val = {}

    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                pos = int(row["A_msa_col"]) +1
                val = float(row["feature_importance_vals_A"])
                interaction = row["interaction"].strip()
            except (ValueError, KeyError):
                continue

            # skip non-interacting rows for normalization/coloring purposes
            if interaction == "non-interacting":
                continue

            # track the max value seen per position across the WHOLE file
            if pos not in full_vals or val > full_vals[pos]:
                full_vals[pos] = val

            # only keep positions matching the requested interaction filter
            if interaction_type == "all" or interaction == interaction_type:
                if pos not in pos_to_val or val > pos_to_val[pos]:
                    pos_to_val[pos] = val

    if not full_vals:
        print("ERROR: No data loaded")
        return

    print(">>> Loaded", len(full_vals), "positions total in dataset")
    print(">>> ", len(pos_to_val), f"positions match interaction='{interaction_type}'")

    # 2. Normalize using the FULL dataset's min/max, so the color scale
    #    is consistent no matter which interaction subset you're viewing
    all_vals = list(full_vals.values())
    vmin, vmax = min(all_vals), max(all_vals)

    # 3. Find every residue number that actually exists in the model,
    #    so anything not in pos_to_val is explicitly zeroed out
    all_positions = get_all_residue_numbers(session, model)
    if not all_positions:
        print(">>> WARNING: could not resolve residues for model", model,
              "- falling back to writing only TSV positions")
        all_positions = set(full_vals.keys())

    # 4. Write attribute file: matched positions get their normalized value,
    #    every other residue in the model gets 0
    attr_file = "/tmp/gini_attribute.defattr"
    with open(attr_file, "w") as out:
        out.write("attribute: gini\nrecipient: residues\n\n")
        for pos in sorted(all_positions):
            if pos in pos_to_val:
                val = pos_to_val[pos]
                norm_val = (val - vmin) / (vmax - vmin + 1e-12)
            else:
                norm_val = 0.0
            out.write(f"\t:{pos}\t{norm_val}\n")

    print(">>> Attribute file written:", attr_file,
          f"({len(all_positions)} residues total, {len(pos_to_val)} colored)")

    # 5. Load attribute and apply coloring
    run(session, f"open {attr_file}")
    run(session, f"color byattribute gini {model} palette white:{palette_color}")
    print(">>> Coloring complete!")


tsv_file = sys.argv[1]
model_id = sys.argv[2]
interaction_type = sys.argv[3] if len(sys.argv) > 3 else "all"  # e.g. "inter_capsomer", "intra_capsomer", "self", or "all"
palette_color = "#a70088"  # change the palette here

run_script(session, tsv_file=tsv_file, model=model_id,
           interaction_type=interaction_type, palette_color=palette_color)