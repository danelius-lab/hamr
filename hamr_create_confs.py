import logging
def main(input_pdbs, output_path, dihedral_atoms, angle_step, prefix, is_amide):
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    count = 0
    for input_pdb in input_pdbs:
        log.info(f"Starting HAMR conformer generation for {input_pdb}.")
        count = generate_conformers(
            input_pdb=input_pdb,
            output_path=output_path,
            dihedral_atoms=dihedral_atoms,
            angle_step=angle_step,
            count=count,
            prefix=prefix,
            log=log,
            is_amide=is_amide
        )
    return

def generate_conformers(input_pdb, output_path, dihedral_atoms, angle_step, count, prefix, log, is_amide):
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    orig_angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *dihedral_atoms)
    iter_start = 0
    iter_end = 361
    if is_amide:
        iter_start = 180 - angle_step
        iter_end = 181 + angle_step
    for angle in range(iter_start, iter_end, angle_step):
        new_angle = orig_angle + angle
        try:
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *dihedral_atoms, new_angle)
            output_file_path = f"{output_path}/{prefix}_{count}.pdb"

            pdb_header = extract_pdb_header(input_pdb=input_pdb, log=log)
            with open(output_file_path, "w") as f:
                pdb_string = pdb_header["full_header"]
                pdb_string += "\n"
                pdb_string += Chem.MolToPDBBlock(mol)
                f.write(pdb_string)
                f.close()
            count += 1
        except:
            log.warn(f"Failed generating conformer for {input_pdb}. Continuing.")
            continue
    return count

def extract_pdb_header(input_pdb, log):
    import re
    crystal_regex = r"CRYST.*\n"
    try:
        with open(input_pdb, "r") as f:
            pdb_string = f.read()
            match = re.findall(crystal_regex, pdb_string)[0]
            split_match = match.split(" ")
            split_match.pop(0)
            unit_cell_dimensions = []
            space_group = []
            while len(unit_cell_dimensions) < 6:
                popped = split_match.pop(0)
                if len(popped) < 1:
                    continue
                unit_cell_dimensions.append(popped.strip())
            unit_cell_dimensions = " ".join(unit_cell_dimensions)
            while len(split_match) > 1:
                popped = split_match.pop(0)
                if len(popped) < 1:
                    continue
                space_group.append(popped.strip())
            space_group = " ".join(space_group)
            f.close()
            return {"full_header": match.strip(), "space_group": space_group, "unit_cell_dimensions": unit_cell_dimensions}

    except Exception as e:
        #TODO: add alternate LLG ranking path instead of exiting
        log.warn(f"Failed extracting PDB header with unit cell dimensions and space group information from this file {input_pdb}. Exiting.")
        exit(1)
        
    
