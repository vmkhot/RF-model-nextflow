from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import csv, sys
import numpy as np


def get_all_residue_numbers(session, model_spec):
    """
    Return a sorted set of all residue numbers (integers) that exist
    in the model matching model_spec, across all chains.
    """
    from chimerax.core.commands import AtomSpecArg, run
    # Resolve the model spec to actual atomic structures
    spec = run(session, f"select {model_spec}", log=False)
    all_positions = set()

    for m in session.models:
        if isinstance(m, AtomicStructure):
            # Only include models that match the requested spec's id
            # (simple check: model id string starts with the requested one)
            if model_spec.lstrip("#") == "" or str(m.id_string) == model_spec.lstrip("#"):
                for r in m.residues:
                    all_positions.add(r.number)

    run(session, "select clear", log=False)
    return all_positions


def run_script(session, tsv_file=None, model="#1", palette_color="bf22a2ff"):
    if tsv_file is None:
        print("ERROR: Must provide tsv_file=/path/to/file")
        return

    print(">>> Loading:", tsv_file)

    # 1. Read feature -> max importance value
    pos_to_val = {}
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                p = int(float(row["feature"].split("_")[0]) + 1)
                v = float(row["feature_importance_vals"])
                if p not in pos_to_val or v > pos_to_val[p]:  # keep max
                    pos_to_val[p] = v
            except:
                pass

    if not pos_to_val:
        print("ERROR: No data loaded")
        return

    print(">>> Loaded", len(pos_to_val), "positions (max value per position)")

    # 2. Normalize
    vals = list(pos_to_val.values())
    vmin, vmax = min(vals), max(vals)

    # 3. Find every residue number that actually exists in the model,
    #    so we can explicitly zero out anything not present in the TSV.
    #    This prevents stale attribute values from a previous run
    #    (e.g. leftover /tmp/gini_attribute.defattr) from sticking around.
    all_positions = get_all_residue_numbers(session, model)
    if not all_positions:
        print(">>> WARNING: could not resolve residues for model", model,
              "- falling back to writing only TSV positions")
        all_positions = set(pos_to_val.keys())

    # 4. Write attribute file: TSV positions get their normalized value,
    #    every other residue in the model gets 0.
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
          f"({len(all_positions)} residues total, {len(pos_to_val)} from TSV)")

    # 5. Load attribute and apply coloring
    run(session, f"open {attr_file}")
    run(session, f"color byattribute gini {model} palette white:{palette_color}")

    print(">>> Coloring complete!")

tsv_file = sys.argv[1]
model_id = sys.argv[2]
palette_color = "#a70088" # change the palette here
run_script(session, tsv_file=tsv_file, model=model_id, palette_color=palette_color)