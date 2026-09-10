# from chimerax.core.commands import run
# import csv, sys
# import numpy as np

# def run_script(session, tsv_file=None, model="#1"):
#     if tsv_file is None:
#         print("ERROR: Must provide tsv_file=/path/to/file")
#         return

#     print(">>> Loading:", tsv_file)

#     # 1. Read seq_position → gini_imp into dict
#     pos_to_val = {}
#     with open(tsv_file) as f:
#         reader = csv.DictReader(f, delimiter="\t")
#         for row in reader:
#             try:
#                 p = int(float(row["seq_position"]))
#                 v = float(row["gini_imp"])
#                 pos_to_val[p] = v
#             except:
#                 pass

#     if not pos_to_val:
#         print("ERROR: No data loaded")
#         return

#     print(">>> Loaded", len(pos_to_val), "positions")

#     # --- Compute min and max
#     vals = list(pos_to_val.values())
#     vmin, vmax = min(vals), max(vals)

#     # --- Convert to ChimeraX attribute file format with normalization
#     attr_file = "/tmp/gini_attribute.defattr"
#     with open(attr_file, "w") as out:
#         out.write("attribute: gini\nrecipient: residues\n\n")
#         for pos, val in pos_to_val.items():
#             norm_val = (val - vmin) / (vmax - vmin + 1e-12)  # normalize 0–1
#             out.write(f"\t:{pos}\t{norm_val}\n")

#     print(">>> Attribute file written:", attr_file)

#     # 3. Load attribute into ChimeraX
#     run(session, f"open {attr_file}")

#     # 4. Apply color mapping
#     run(session, f"color byattribute gini {model} palette Oranges")

#     print(">>> Coloring complete!")
    
# tsv_file = sys.argv[1]
# model_id = sys.argv[2]
# run_script(session, tsv_file=tsv_file, model=model_id)

# from chimerax.core.commands import run
# import csv, sys

# def run_script(session, tsv_file=None, model="#1"):
#     if tsv_file is None:
#         print("ERROR: Must provide tsv_file=/path/to/file")
#         return

#     print(">>> Loading:", tsv_file)

#     # 1. Read seq_position → max gini_imp
#     pos_to_val = {}
#     with open(tsv_file) as f:
#         reader = csv.DictReader(f, delimiter="\t")
#         for row in reader:
#             try:
#                 # p = int(float(row["seq_position"]))
#                 # v = float(row["gini_imp"])
#                 p = int(float(row["feature"].split("_")[0])+1)
#                 v = float(row["feature_importance_vals"])
#                 if p not in pos_to_val or v > pos_to_val[p]:  # keep max
#                     pos_to_val[p] = v
#             except:
#                 pass

#     if not pos_to_val:
#         print("ERROR: No data loaded")
#         return

#     print(">>> Loaded", len(pos_to_val), "positions (max gini per position)")

#     # 2. Normalize and write attribute file
#     vals = list(pos_to_val.values())
#     vmin, vmax = min(vals), max(vals)

#     attr_file = "/tmp/gini_attribute.defattr"
#     with open(attr_file, "w") as out:
#         out.write("attribute: gini\nrecipient: residues\n\n")
#         for pos, val in pos_to_val.items():
#             norm_val = (val - vmin) / (vmax - vmin + 1e-12)
#             out.write(f"\t:{pos}\t{norm_val}\n")

#     print(">>> Attribute file written:", attr_file)

#     # 3. Load attribute and apply coloring
#     run(session, f"open {attr_file}")
#     run(session, f"color byattribute gini {model} palette Oranges")

#     print(">>> Coloring complete!")

# tsv_file = sys.argv[1]
# model_id = sys.argv[2]
# run_script(session, tsv_file=tsv_file, model=model_id)

from chimerax.core.commands import run
import csv, sys
import numpy as np


def run_script(session, tsv_file=None, model="#1", palette_color="bf22a2ff"):
    if tsv_file is None:
        print("ERROR: Must provide tsv_file=/path/to/file")
        return

    print(">>> Loading:", tsv_file)

    # 1. Read feature → max importance value
    pos_to_val = {}
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # if row['biome_predictor'] == "ocean":
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

    # 2. Normalize and write attribute file
    vals = list(pos_to_val.values())
    vmin, vmax = min(vals), max(vals)

    attr_file = "/tmp/gini_attribute.defattr"
    with open(attr_file, "w") as out:
        out.write("attribute: gini\nrecipient: residues\n\n")
        for pos, val in pos_to_val.items():
            norm_val = (val - vmin) / (vmax - vmin + 1e-12)
            # norm_val = val / vmax
            out.write(f"\t:{pos}\t{norm_val}\n")

    print(">>> Attribute file written:", attr_file)

    # 3. Load attribute and apply coloring
    try:
        run(session, "delattr gini residues")
    except Exception:
        pass
    run(session, f"open {attr_file}")
    run(session, f"color byattribute gini {model} palette white:{palette_color}")

    print(">>> Coloring complete!")

tsv_file = sys.argv[1]
model_id = sys.argv[2]
palette_color = "#a70088" # change the palette here
run_script(session, tsv_file=tsv_file, model=model_id, palette_color=palette_color)