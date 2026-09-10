from pymol import cmd, util
import pandas as pd
import matplotlib
import matplotlib.cm as cm

def color_by_gini(table_file, obj_name="protein", cmap_name="Greens"):
    """
    Color residues in PyMOL based on Gini importance values.

    Parameters
    ----------
    table_file : str
        Path to the input TSV or CSV file containing columns:
        ['index', 'shap_value', 'gini_imp', 'cluster_num', 'msa_position', 'AA', 'seq_position']
    obj_name : str
        The name of the PyMOL object (already loaded) to color.
    cmap_name : str
        Matplotlib colormap name (default: 'bwr' for blue-white-red)
    """

    # --- Load table
    df = pd.read_csv(table_file, sep=None, engine="python")  # auto-detects comma or tab
    if "seq_position" not in df.columns or "gini_imp" not in df.columns:
        print("Error: Table must contain 'seq_position' and 'gini_imp' columns.")
        return

    # --- Prepare data
    df = df.dropna(subset=["seq_position", "gini_imp"])
    df = df.drop_duplicates("seq_position")
    df["seq_position"] = df["seq_position"].astype(int)

    # Normalize gini values for coloring
    gini_min, gini_max = df["gini_imp"].min(), df["gini_imp"].max()
    df["norm"] = (df["gini_imp"] - gini_min) / (gini_max - gini_min + 1e-12)
    print(df[["seq_position", "gini_imp", "norm"]].head())

    cmap = cm.get_cmap(cmap_name)

    # --- Create colors in PyMOL
    for _, row in df.iterrows():
        resi = int(row["seq_position"])
        color_rgb = cmap(row["norm"])[:3]
        color_name = f"gini_{resi}"
        cmd.set_color(color_name, color_rgb)
        cmd.color(color_name, f"{obj_name} and resi {resi}")

    # Optional: pretty-up the view
    # util.cbc(obj_name)
    cmd.show("cartoon", obj_name)
    cmd.bg_color("black")
    cmd.orient(obj_name)
    print(f"✅ Colored {len(df)} residues in {obj_name} using '{cmap_name}' colormap.")

# Make the function available in PyMOL
cmd.extend("color_by_gini", color_by_gini)