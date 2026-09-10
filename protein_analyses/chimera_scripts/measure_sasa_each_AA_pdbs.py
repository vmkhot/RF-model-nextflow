from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import os

# ===== USER INPUT =====
monomer_dir = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/Microviridae_analysis/all_proteins_rank001_folded/rank001"
output_file = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/mapping_monomers_BSA/sasa_per_residue.tsv"

# One-letter code lookup
THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O','ASX':'B','GLX':'Z','UNK':'X'
}

def get_atomic_models(session):
    return [m for m in session.models.list() if isinstance(m, AtomicStructure)]

def get_per_residue_sasa(session, model):
    """
    Run measure sasa, then read per-atom area attribute back
    and sum by residue.
    """
    run(session, f"measure sasa #{model.id_string}")

    residue_sasa = {}
    for atom in model.atoms:
        res = atom.residue
        key = (res.number, res.name)
        area = getattr(atom, 'area', 0.0) or 0.0
        residue_sasa[key] = residue_sasa.get(key, 0.0) + area

    return residue_sasa

def process_file(session, pdb_path):
    models_before = set(m.id_string for m in get_atomic_models(session))

    run(session, f"open {pdb_path}")

    new_models = [m for m in get_atomic_models(session)
                  if m.id_string not in models_before]
    if not new_models:
        raise RuntimeError("No model loaded")

    model = new_models[0]
    residue_sasa = get_per_residue_sasa(session, model)
    run(session, f"close #{model.id_string}")

    return residue_sasa

def run_all(session):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    pdb_files = sorted([f for f in os.listdir(monomer_dir) if f.endswith(".pdb")])
    print(f"Found {len(pdb_files)} PDB files\n")

    with open(output_file, "w") as out:
        out.write("protein\tposition\tAA\tsasa\n")

        for i, file in enumerate(pdb_files, 1):
            print(f"[{i}/{len(pdb_files)}] {file}")
            path = os.path.join(monomer_dir, file)
            protein_name = os.path.splitext(file)[0]

            try:
                residue_sasa = process_file(session, path)

                # Sort by residue number and write rows
                for (res_num, res_name), sasa in sorted(residue_sasa.items()):
                    one_letter = THREE_TO_ONE.get(res_name.upper(), 'X')
                    out.write(f"{protein_name}\t{res_num}\t{one_letter}\t{sasa:.4f}\n")

                print(f"    Written {len(residue_sasa)} residues")

            except Exception as e:
                print(f"    ERROR: {e}")
                out.write(f"{protein_name}\tERROR\tERROR\tERROR\n")

            out.flush()

    print("\nAll done!")

# ===== ENTRY POINT =====
if 'session' in globals():
    run_all(session)