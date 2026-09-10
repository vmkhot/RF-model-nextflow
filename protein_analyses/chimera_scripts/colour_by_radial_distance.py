# this script calculates the distance from the center of the capsid to the outer surface and saves as an attribute (radial_dist) and colours is blue to red... this is to see what is on the inner surface and what is on the outer most


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

run(session, "color byattribute radial_dist #2 palette blue:white:red")