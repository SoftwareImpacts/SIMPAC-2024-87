import json
import shutil
import sys
from core_functions import *
from pysimm_system import PysimmSystem

if __name__ == "__main__":
    log_message_to_file(message=f"Program initialized.")
    output_dir = os.path.join(settings.output_dir, "")
    os.makedirs(output_dir)
    pdb_file_path = settings.pdb_file_path
    pdb_file_lines = read_pdb_file(pdb_file_path)
    log_message_to_file(message=f"{pdb_file_path} file read.")

    box_dimensions = extract_pdb_box_dimensions(pdb_header=pdb_file_lines[0])
    pdb_molecule_blocks = extract_molecules(pdb_file_lines)
    log_message_to_file(message=f"{len(pdb_molecule_blocks)} molecules detected.")

    smiles_all = smiles_from_blocks(pdb_molecule_blocks)
    log_message_to_file(message=f"SMILES notations generated successfully.")

    molecule_counts, molecule_objects_by_mol_type = aggregate_molecules(smiles_all)

    mol_types_map = smiles_to_mol_type_map(dictionary=molecule_counts)
    molecular_properties = compute_molecular_properties(molecule_counts)
    log_message_to_file(
        message=f"Molecular properties computed: {json.dumps(molecular_properties, indent=4)}"
    )

    mix_properties = mixture_properties(molecular_properties)
    log_message_to_file(
        message=f"Mixture properties computed: {json.dumps(mix_properties, indent=4)}"
    )

    molecule_objects_by_mol_type_nbr = convert_smiles_type_to_number(
        dictionary=molecule_objects_by_mol_type, map_dictionary=mol_types_map
    )

    pysimm_system = PysimmSystem(
        molecule_objects_by_mol_type_nbr=molecule_objects_by_mol_type_nbr,
        force_field=settings.force_field,
        charges=settings.charges,
    )

    pysimm_system.initialize_system()
    log_message_to_file(message=f"Pysimm system created.")

    pysimm_system.load_all_molecules_into_system()
    log_message_to_file(
        message=f"Force field [{settings.force_field}] types and charges [{settings.charges}] assigned."
    )
    pysimm_system.set_box_dimensions(dimensions=box_dimensions)
    log_message_to_file(
        message=f"PDB dimensions [{json.dumps(box_dimensions, indent=4)}] assigned to Pysimm system."
    )

    pysimm_system.generate_lammps_inputs()
    log_message_to_file(message=f"LAMMPS structure.dat file generated.")

    log_message_to_file(message=f"Program terminated successfully.")
    shutil.move(src="calls_log.txt", dst=os.path.join(output_dir, f"calls_log.txt"))
    os.remove(path="pysimm.sim.in")
    sys.exit()
