from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import os

# ===== USER INPUT =====
monomer_dir   = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/Microviridae_analysis/all_proteins_rank001_folded/rank001"
output_file   = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/mapping_monomers_BSA/sasa_results.tsv"


def get_atomic_models(session):
    return [m for m in session.models.list() if isinstance(m, AtomicStructure)]


def get_sasa(session, model):
    """Robust SASA extraction across ChimeraX versions"""
    import re
    import io
    from contextlib import redirect_stdout

    # Try direct return
    result = run(session, f"measure sasa #{model.id_string}")
    if isinstance(result, (int, float)):
        return float(result)

    # Try atom attributes
    for attr in ["ses_areas", "accessibilities", "sasa_areas", "solvent_areas"]:
        try:
            return float(getattr(model.atoms, attr).sum())
        except Exception:
            pass

    # Fallback: parse log output
    buf = io.StringIO()
    with redirect_stdout(buf):
        run(session, f"measure sasa #{model.id_string}")
    output = buf.getvalue()

    m = re.search(r'=\s*([\d.]+)', output)
    if m:
        return float(m.group(1))

    raise RuntimeError(f"Could not retrieve SASA for #{model.id_string}")


def process_file(session, pdb_path):
    models_before = set(m.id_string for m in get_atomic_models(session))

    # Load structure
    run(session, f"open {pdb_path}")

    new_models = [m for m in get_atomic_models(session)
                  if m.id_string not in models_before]
    if not new_models:
        raise RuntimeError("No model loaded")

    model = new_models[0]

    # Compute SASA
    sasa = get_sasa(session, model)

    # Clean up
    run(session, f"close #{model.id_string}")

    return sasa


def run_all(session):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    pdb_files = sorted([f for f in os.listdir(monomer_dir) if f.endswith(".pdb")])
    print(f"Found {len(pdb_files)} PDB files\n")

    with open(output_file, "w") as out:
        out.write("protein\tsasa\n")

        for i, file in enumerate(pdb_files, 1):
            print(f"[{i}/{len(pdb_files)}] {file}")
            path = os.path.join(monomer_dir, file)

            try:
                sasa = process_file(session, path)
                print(f"    SASA = {sasa:.2f}")
                out.write(f"{file}\t{sasa:.4f}\n")
            except Exception as e:
                print(f"    ERROR: {e}")
                out.write(f"{file}\tERROR\n")

            out.flush()

    print("\nAll done!")


# ===== ENTRY POINT =====
if 'session' in globals():
    run_all(session)