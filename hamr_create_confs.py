import logging
from sys import exit

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
    from rdkit.Chem import rdMolTransforms, rdmolops

    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    orig_angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *dihedral_atoms["angle_indices"])
    iter_start = 0
    iter_end = 361
    if not is_order_correct(input_pdb=input_pdb, atom_idx1=dihedral_atoms['angle_indices'][1], atom_idx2=dihedral_atoms['angle_indices'][2], log=log):
            dihedral_atoms['angle_indices'] = dihedral_atoms['angle_indices'][::-1]
            dihedral_atoms['angle_atom_names'] = dihedral_atoms['angle_atom_names'][::-1]
    if is_amide:
        iter_start = 170
        iter_end = 191

        
        atom_idx_dict = {}
        for atom in mol.GetAtoms():
                atom_idx_dict[atom.GetMonomerInfo().GetName().strip()] = atom.GetIdx()
                atom_idx_dict[str(atom.GetIdx())] = atom.GetMonomerInfo().GetName().strip()

        for neighbor in mol.GetAtomWithIdx(atom_idx_dict[dihedral_atoms["angle_atom_names"][1]]).GetNeighbors():
                if "C" in dihedral_atoms["angle_atom_names"][1]:
                    if neighbor.GetAtomicNum() == 8:
                        dihedral_atoms["angle_atom_names"][0] = atom_idx_dict[str(neighbor.GetIdx())]
                        dihedral_atoms["angle_indices"][0] = neighbor.GetIdx()
                elif "N" in dihedral_atoms["angle_atom_names"][1]:
                    if neighbor.GetAtomicNum() == 1:
                       dihedral_atoms["angle_atom_names"][0] = atom_idx_dict[str(neighbor.GetIdx())]
                       dihedral_atoms["angle_indices"][0] = neighbor.GetIdx()
                        
        for neighbor in mol.GetAtomWithIdx(atom_idx_dict[dihedral_atoms["angle_atom_names"][2]]).GetNeighbors():
                if "C" in dihedral_atoms["angle_atom_names"][2]:
                    if neighbor.GetAtomicNum() == 8:
                        dihedral_atoms["angle_atom_names"][3] = atom_idx_dict[str(neighbor.GetIdx())]
                        dihedral_atoms["angle_indices"][3] = neighbor.GetIdx()
                elif "N" in dihedral_atoms["angle_atom_names"][2]:
                    if neighbor.GetAtomicNum() == 1:
                        dihedral_atoms["angle_atom_names"][3] = atom_idx_dict[str(neighbor.GetIdx())]
                        dihedral_atoms["angle_indices"][3] = neighbor.GetIdx()
    for angle in range(iter_start, iter_end, angle_step):
        new_angle = orig_angle + angle
        try:
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *dihedral_atoms["angle_indices"], new_angle)
            output_file_path = f"{output_path}/{prefix}_{count}.pdb"
            if is_amide:
                dist_matrix = rdmolops.Get3DDistanceMatrix(mol)
                print("AMIDE ATOM DISTANCE")
                print(dist_matrix[dihedral_atoms['angle_indices'][0]][dihedral_atoms['angle_indices'][-1]])
                if dist_matrix[dihedral_atoms['angle_indices'][0]][dihedral_atoms['angle_indices'][-1]] < 3:
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(), *dihedral_atoms["angle_indices"], 180 + new_angle)
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
        
    
def is_order_correct(input_pdb, atom_idx1, atom_idx2, log):
    rank = morgan_algo(input_pdb=input_pdb, num_cycles=3, log=log)
    return rank[atom_idx1] > rank[atom_idx2]


def morgan_algo(input_pdb, num_cycles, log):
    from rdkit import Chem
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    rank = []
    for atom in mol.GetAtoms():
        rank.append(1)
    while num_cycles > 0:
        new_rank = [*rank]
        for atom in mol.GetAtoms():
            new_rank_entry = new_rank[atom.GetIdx()]
            for neighbor in atom.GetNeighbors():
                new_rank_entry += rank[neighbor.GetIdx()]
            new_rank[atom.GetIdx()] = new_rank_entry
        num_cycles -= 1
        rank = new_rank
    return rank
