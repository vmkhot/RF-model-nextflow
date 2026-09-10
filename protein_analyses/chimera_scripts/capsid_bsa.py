from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import os

# ===== USER INPUT =====
template_path = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/8DES-assembly1.cif"
monomer_dir   = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/Microviridae_analysis/multiple_cluster_analysis/mixed_cluster_modelling/example_pbds/pdbs"
output_file   = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/mapping_monomers_BSA/bsa_results.tsv"


def get_atomic_models(session):
    """Return only AtomicStructure models, ignoring PseudobondGroups etc."""
    return [m for m in session.models.list() if isinstance(m, AtomicStructure)]


def get_sasa(session, model):
    """Compute SASA for a model - compatible with ChimeraX 1.10.x"""
    import re
    
    # measure sasa returns the area value directly in some versions
    result = run(session, f"measure sasa #{model.id_string}")
    
    # Try result directly (float in some versions)
    if isinstance(result, (int, float)):
        return float(result)
    
    # Try atoms attributes (attribute name varies by version)
    for attr in ["ses_areas", "accessibilities", "sasa_areas", "solvent_areas"]:
        try:
            return float(getattr(model.atoms, attr).sum())
        except AttributeError:
            continue

    # Last resort: parse from the printed log line we can already see working:
    # "Solvent accessible area for #11 (4249 atoms) = 34307"
    import io
    from contextlib import redirect_stdout
    buf = io.StringIO()
    with redirect_stdout(buf):
        run(session, f"measure sasa #{model.id_string}")
    output = buf.getvalue()
    m = re.search(r'Solvent accessible area.*?=\s*([\d.]+)', output)
    if m:
        return float(m.group(1))

    raise RuntimeError(f"Could not retrieve SASA for #{model.id_string}")


def process_monomer(session, monomer_path, template):
    models_before = set(m.id_string for m in get_atomic_models(session))

    # 1. Load monomer
    run(session, f"open {monomer_path}")
    new_models = [m for m in get_atomic_models(session)
                  if m.id_string not in models_before]
    if not new_models:
        raise RuntimeError(f"No new atomic model loaded for {monomer_path}")
    monomer = new_models[0]
    print(f"    Monomer loaded as #{monomer.id_string} ({len(list(monomer.chains))} chains)")

    # 2. Align monomer to first chain of template
    first_chain = list(template.chains)[0]
    run(session, f"matchmaker #{monomer.id_string} to #{template.id_string}/{first_chain.chain_id}")

    # 3. Compute isolated monomer SASA
    sasa_monomer = get_sasa(session, monomer)

    # 4. Build full capsid - only use 'A' chains from template (major capsid protein)
    # a_chains = [c for c in template.chains if c.chain_id == "A"]
    a_chains = [c for c in template.chains if c.structure is template and c.polymer_type is not None]
    print(f"    Building capsid from {len(a_chains)} A-chains...")

    capsid_models = []
    for chain in a_chains:
        models_before_copy = set(m.id_string for m in get_atomic_models(session))
        run(session, f"combine #{monomer.id_string} close false")
        copy_model = [m for m in get_atomic_models(session)
                      if m.id_string not in models_before_copy][0]
        run(session, f"matchmaker #{copy_model.id_string} to #{template.id_string}/{chain.chain_id}")
        capsid_models.append(copy_model)

    # 5. Merge all copies into one capsid model
    model_ids = " ".join([f"#{m.id_string}" for m in capsid_models])
    models_before_capsid = set(m.id_string for m in get_atomic_models(session))
    run(session, f"combine {model_ids} name capsid_model close false")
    capsid = [m for m in get_atomic_models(session)
              if m.id_string not in models_before_capsid][0]

    # 6. Compute capsid SASA
    sasa_capsid = get_sasa(session, capsid)

    # 7. BSA
    n   = len(capsid_models)
    bsa = (n * sasa_monomer) - sasa_capsid

    # 8. Clean up
    run(session, f"close #{monomer.id_string}")
    for m in capsid_models:
        try:
            run(session, f"close #{m.id_string}")
        except Exception:
            pass
    run(session, f"close #{capsid.id_string}")

    return sasa_monomer, sasa_capsid, bsa, n


def run_all(session):
    print(">>> Loading template...")
    run(session, f"open {template_path}")

    # Get only atomic structures — avoids grabbing PseudobondGroup
    atomic_models = get_atomic_models(session)
    if not atomic_models:
        print("ERROR: No atomic models found after loading template")
        return

    template = atomic_models[0]
    all_chains = list(template.chains)
    a_chains   = [c for c in all_chains if c.chain_id.startswith("A")]
    print(f">>> Template loaded: {template.name}")
    print(f">>> Total chains: {len(all_chains)}, A-chains (MCP): {len(a_chains)}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    pdb_files = sorted([f for f in os.listdir(monomer_dir) if f.endswith(".pdb")])
    print(f">>> Found {len(pdb_files)} PDB files to process\n")

    with open(output_file, "w") as out:
        out.write("protein\tsasa_monomer\tsasa_capsid\tbsa\tn_chains\n")

        for i, file in enumerate(pdb_files, 1):
            print(f">>> [{i}/{len(pdb_files)}] Processing {file}...")
            path = os.path.join(monomer_dir, file)
            try:
                sasa_m, sasa_c, bsa, n = process_monomer(session, path, template)
                print(f"    SASA monomer={sasa_m:.2f}  capsid={sasa_c:.2f}  BSA={bsa:.2f}  n={n}")
                out.write(f"{file}\t{sasa_m:.4f}\t{sasa_c:.4f}\t{bsa:.4f}\t{n}\n")
                out.flush()
            except Exception as e:
                print(f"    ERROR: {e}")
                out.write(f"{file}\tERROR\tERROR\tERROR\tERROR\n")
                out.flush()

    print("\n>>> All done!")


run_all(session)