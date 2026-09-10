from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import os

# ===== USER INPUT =====
monomer_mmcif = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/8des_gokushovirus_ec6098.cif"
capsid_mmcif  = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/8DES-assembly1.cif"
output_dir    = "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/reference_pdbs/mapping_monomers_BSA/"

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O','ASX':'B','GLX':'Z','UNK':'X'
}

def get_atomic_models(session):
    return [m for m in session.models.list() if isinstance(m, AtomicStructure)]

def get_per_residue_sasa(session, model_spec):
    """Run SASA on a model spec and return {(res_num, chain_id): (aa_one, sasa)}"""
    run(session, f"measure sasa {model_spec}")

    # Retrieve the model object to iterate atoms
    all_models = get_atomic_models(session)

    # Find the model matching the spec (by id_string prefix)
    target_id = model_spec.lstrip("#")
    residue_sasa = {}

    for model in all_models:
        if not model.id_string.startswith(target_id):
            continue
        for atom in model.atoms:
            res      = atom.residue
            key      = (res.number, res.chain_id)
            area     = getattr(atom, 'area', 0.0) or 0.0
            if key not in residue_sasa:
                residue_sasa[key] = {'aa': THREE_TO_ONE.get(res.name.upper(), 'X'), 'sasa': 0.0}
            residue_sasa[key]['sasa'] += area

    return residue_sasa

def write_sasa(residue_sasa, output_path, protein_name, chain_filter=None):
    """Write per-residue SASA to TSV, optionally filtering to one chain."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("protein\tposition\tchain\tAA\tsasa\n")
        for (res_num, chain_id), vals in sorted(residue_sasa.items()):
            if chain_filter and chain_id != chain_filter:
                continue
            f.write(f"{protein_name}\t{res_num}\t{chain_id}\t{vals['aa']}\t{vals['sasa']:.4f}\n")
    print(f"  Written: {output_path}")

def run_all(session):
    os.makedirs(output_dir, exist_ok=True)

    # ── 1. Monomer ────────────────────────────────────────────────────────────
    print("Loading monomer...")
    models_before = set(m.id_string for m in get_atomic_models(session))
    run(session, f"open {monomer_mmcif}")
    new_models = [m for m in get_atomic_models(session)
                  if m.id_string not in models_before]
    if not new_models:
        raise RuntimeError("Monomer failed to load")

    monomer_model = new_models[0]
    monomer_spec  = f"#{monomer_model.id_string}"
    print(f"  Loaded as {monomer_spec}")

    monomer_sasa = get_per_residue_sasa(session, monomer_spec)
    write_sasa(
        monomer_sasa,
        os.path.join(output_dir, "sasa_per_residue_monomer.tsv"),
        protein_name  = "EC6098_monomer",
        chain_filter  = None   # monomer has one chain; keep all
    )
    run(session, f"close {monomer_spec}")

    # ── 2. Capsid assembly — chain A-53 only ─────────────────────────────────
    print("\nLoading capsid assembly...")
    models_before = set(m.id_string for m in get_atomic_models(session))
    run(session, f"open {capsid_mmcif}")
    new_models = [m for m in get_atomic_models(session)
                  if m.id_string not in models_before]
    if not new_models:
        raise RuntimeError("Capsid failed to load")

    capsid_model = new_models[0]
    capsid_spec  = f"#{capsid_model.id_string}"
    print(f"  Loaded as {capsid_spec}")

    # Run SASA on the FULL assembly (so buried surface is computed correctly
    # in the context of all neighbouring chains), then extract chain A-53 only
    print("  Computing SASA on full assembly (this may take a while)...")
    capsid_sasa = get_per_residue_sasa(session, capsid_spec)

    # Find the chain ID that corresponds to A-53 in the mmcif
    # ChimeraX names chains from mmcif auth_asym_id; A-53 is typically "A-53"
    # Adjust chain_filter below if ChimeraX renames it (check Log after opening)
    target_chain = "A-53"
    write_sasa(
        capsid_sasa,
        os.path.join(output_dir, "sasa_per_residue_capsid.tsv"),
        protein_name  = "EC6098_capsid_A53",
        chain_filter  = target_chain
    )
    run(session, f"close {capsid_spec}")

    print("\nAll done!")

# ===== ENTRY POINT =====
if 'session' in globals():
    run_all(session)