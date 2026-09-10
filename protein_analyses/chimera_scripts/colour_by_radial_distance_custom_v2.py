# this script calculates the distance from the center of the capsid to the outer surface and saves as an attribute (radial_dist) and colours it in 3 discrete bins:
#   < 115 A        : #4575b4ff (blue)
#   115 - 120 A    : #ffc117ff (yellow)
#   >= 120 A       : #d73027ff (red)

import numpy as np
from chimerax.atomic import AtomicStructure, Residue
from chimerax.core.commands import run

# Get the atomic structure specifically, not surface objects
structure = [m for m in session.models if isinstance(m, AtomicStructure) and m.id_string == "2"][0]

coords = structure.atoms.coords
center = coords.mean(axis=0)
print("Center:", center)

for res in structure.residues:
    ca = [a for a in res.atoms if a.name == "CA"]
    if ca:
        res.radial_dist = float(np.linalg.norm(ca[0].coord - center))

Residue.register_attr(session, "radial_dist", "radial coloring", attr_type=float)

# Build a palette with control points placed right at the bin edges (with a
# tiny epsilon gap) so the color jumps discretely instead of blending smoothly.
# ChimeraX clamps values outside the given range to the nearest endpoint color,
# so anything < 115 gets blue and anything >= 120 gets red automatically.
# Boundaries placed so that:
#   value <  115        -> blue
#   115 <= value < 120  -> yellow
#   value >= 120        -> red
eps = 0.001
palette = (
    f"{115-eps},#1B5E20:"
    f"115,#ccc19f:"
    f"{120-eps},#ccc19f:"
    f"120,#7A1F3D"
)

run(session, f"color byattribute radial_dist #2 palette {palette}")

# diverging_cols <- c(
#   low  = "#1B5E20",  # forest green
#   mid  = "#F8F4E8",  # cream
#   high = "#7A1F3D"   # burgundy
# )