import matplotlib.cm as cm
import numpy as np

cmap = cm.get_cmap("Greens")

with open("greens_matplotlib.cmap", "w") as f:
    for i in range(256):
        r, g, b, _ = cmap(i/255)
        f.write(f"{r} {g} {b}\n")