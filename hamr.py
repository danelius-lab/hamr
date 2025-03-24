import logging
import os
import re 
from math import floor

import hamr_create_confs
import hamr_calc_struct_fact
import phase
import refine
import analyze_phase_results
def main(
        input_pdb, 
        input_mtz, 
        restraint_cif, 
        output_path, 
        smiles_string,
        r_free_fraction,
        refinement_columns,
        r_factor_columns,
        num_refine_cycles,
        num_to_persist,
        prefix,
        angle_step,
        #TODO: implement this option (idk if this is actually needed/useful)
        faulty_conformer_rmsd_cutoff,
        should_force_trans_amides,
        fraction_to_phase
        ):
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    amide_bonds = find_amides(input_pdb=input_pdb, template_smiles_string=smiles_string, log=log)
    dihedral_path = rank_dihedrals(input_pdb=input_pdb, angle_step=angle_step, log=log)
    input_pdbs = [{"file_path": input_pdb, "llg": 0, "r_factor": 1, "name": "initial_pdb"}]
    
    #TODO: add time handling/counting messages
    round_count = 0
    for dihedral_angle in dihedral_path:
        is_amide = False
        for amide in amide_bonds:
            print(dihedral_angle)
            print(amide)
            if dihedral_angle["angle_atom_names"][1] in amide and dihedral_angle["angle_atom_names"][2] in amide and should_force_trans_amides:
                is_amide = True
                #TODO: verify this is working as intended
                log.info("Found dihedral angle for amide bond in this HAMR cycle. Forcing trans-planar geometry and continuing.")
                break
        
        round_output_path = f"{output_path}/ROUND_{round_count}"
        os.makedirs(round_output_path, exist_ok=True)

        input_pdb_paths = [x["file_path"] for x in input_pdbs]
        hamr_create_confs.main(
            input_pdbs=input_pdb_paths, 
            output_path=round_output_path, 
            dihedral_atoms=dihedral_angle,
            angle_step=angle_step,
            prefix=prefix,
            is_amide=is_amide)
        
        struct_facts = []
       
        for file in os.listdir(round_output_path):
            if ".pdb" not in file:
                continue
            full_file_path = f"{round_output_path}/{file}"
            struct_fact = hamr_calc_struct_fact.main(
                input_mtz=input_mtz,
                input_pdb=full_file_path,
                intensity_column=r_factor_columns[0],
                sigma_column=r_factor_columns[1],
            )
            conf_num = int(file.split("_")[-1].split(".pdb")[0])
            input_pdb_idx = floor(conf_num / (370/angle_step))
            print(input_pdb_idx)
            orig_struct_fact = hamr_calc_struct_fact.main(
                input_mtz=input_mtz,
                input_pdb=input_pdbs[input_pdb_idx]['file_path'],
                intensity_column=r_factor_columns[0],
                sigma_column=r_factor_columns[1],
            )
            struct_facts.append({"file_path": full_file_path, "struct_fact": struct_fact - orig_struct_fact})

        if round_count == 0:
            input_pdbs = []
        struct_facts = list(sorted(struct_facts, key=lambda x : x["struct_fact"]))
        print(struct_facts)
        phaser_path = f"{round_output_path}/PHASER"
        os.makedirs(phaser_path, exist_ok=True)
        count = 0
        valid_structs = 0
        while (valid_structs < fraction_to_phase * len(struct_facts) or valid_structs < num_to_persist) and count < len(struct_facts):
            prefix_conf = struct_facts[count]["file_path"].split("/")[-1].split(".pdb")[0]
            if not validate_pdb(struct_facts[count]["file_path"], log):
                count += 1
                continue
            if  not should_append_to_list(struct_facts[count]["file_path"], input_pdbs):
                log.info(f"Duplicate conformer already in pointer list detected: {prefix_conf}. Skipping MR for this conformer and continuing.")
                count += 1
                continue
            prefix_conf = struct_facts[count]["file_path"].split("/")[-1].split(".pdb")[0]
            phased_pdb = phase.run_phaser(
                input_mtz=input_mtz,
                input_pdb=struct_facts[count]["file_path"],
                output_path=phaser_path,
                prefix=prefix_conf,
                log=log,
            )
            if validate_pdb(input_pdb=phased_pdb, log=log):
                valid_structs += 1
                solu_file = re.sub(".1.pdb", ".sol", phased_pdb)
                with open(solu_file,"r") as f:
                    stats = analyze_phase_results.extract_solu_stats(solu_string=f.read(), log=log, file_path=phased_pdb, model_name=phased_pdb.split("/")[-1], input_mtz=input_mtz, intensity_column=r_factor_columns[0], sigma_column=r_factor_columns[1])
                input_pdbs.append(stats)
                log.info("Current solutions:")
                for f in input_pdbs:
                    log.info(f"{f['name']} -- LLG: {f['llg']} R-factor: {f['r_factor']}")
            count += 1
        # for phased_file in os.listdir(phaser_path):
        #     if ".pdb" not in phased_file:
        #         continue
        #     input_pdbs.append(f"{phaser_path}/{phased_file}")
        input_pdbs = list(sorted(input_pdbs, key=lambda x: x["llg"], reverse=True))
        input_pdbs = input_pdbs[:num_to_persist]
        round_count += 1
        log.info(f"Finished HAMR cycle {round_count} of {len(dihedral_path)}. Continuing with next cycle.")
    log.info(f"Starting initial refinement of top {num_to_persist} candidates.")
    refine_output_path = f"{output_path}/ROUND_{round_count}"
    os.makedirs(refine_output_path, exist_ok=True)
    return refine.main(
        input_pdbs=input_pdbs,
        input_mtz=input_mtz,
        restraint_cif=restraint_cif,
        num_refine_cycles=num_refine_cycles,
        r_free_fraction=r_free_fraction,
        refinement_columns=refinement_columns,
        output_dir=refine_output_path
    )

def validate_pdb(input_pdb, log):
    from rdkit import Chem
    try:
        test_mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
        if test_mol == None:
            raise Exception()
        return True
    except:
        log.warn(f"Invalid PDB found for {input_pdb}. This structure will be discarded. Continuing.")
        return False

def find_amides(input_pdb, template_smiles_string, log):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    log.info("Starting process to find amides.")
    amide_smarts = "C(=O)-N"
    amide_query = Chem.MolFromSmarts(amide_smarts)

    mol = Chem.MolFromPDBFile(input_pdb)
    template_mol = Chem.MolFromSmiles(template_smiles_string)
    mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
    matches = mol.GetSubstructMatches(amide_query)
    num_matches = len(matches)
    log.info(f"Found {num_matches} amide bonds.")
    
    return [[transform_idx_to_name(match[0], mol),transform_idx_to_name(match[1], mol),transform_idx_to_name(match[2], mol)] for match in matches]



def rank_dihedrals(input_pdb, angle_step, log):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign, rdMolTransforms

    log.info("Starting process to rank all modifiable dihedral angles.")
    dihed_smarts_query = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
    dihed_smarts_mol = Chem.MolFromSmarts(dihed_smarts_query)
    
    mol = Chem.MolFromPDBFile(input_pdb)
    matches = mol.GetSubstructMatches(dihed_smarts_mol)

    dihed_list = []
    for match in matches:
        idx2 = match[0]
        idx3 = match[1]
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        jAtom = mol.GetAtomWithIdx(idx2)
        kAtom = mol.GetAtomWithIdx(idx3)
        for b1 in jAtom.GetBonds():
            if (b1.GetIdx() == bond.GetIdx()):
                continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in kAtom.GetBonds():
                if ((b2.GetIdx() == bond.GetIdx())
                or (b2.GetIdx() == b1.GetIdx())):
                    continue
                idx4 = b2.GetOtherAtomIdx(idx3)
            # skip 3-membered rings
                if (idx4 == idx1):
                    continue
                dihed_list.append([idx1, idx2, idx3, idx4])
    
    #double checking that all dihedral angles are manipulatable and won't raise
    allowed_dihed_angles = []
    for angle in dihed_list:
        try:
            orig_angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3])
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(),angle[0],angle[1],angle[2],angle[3], orig_angle)
        except:
            continue
        allowed_dihed_angles.append(angle)

    ranked_allowed_dihed_angles = []
    ref_mol = Chem.MolFromPDBFile(input_pdb)
    for dihed_angle in allowed_dihed_angles:
        del mol
        mol = Chem.MolFromPDBFile(input_pdb)
        rmsd_list = []
        for angle in range(angle_step, angle_step+360, angle_step):
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(), dihed_angle[0], dihed_angle[1], dihed_angle[2], dihed_angle[3], angle)
            rmsd = rdMolAlign.GetBestRMS(mol, ref_mol)
            if rmsd > 1000:
                log.warn("Invalid dihedral angle found when ranking dihedrals, this angle is most likely colinear, applying penalty to this dihedral angle and skipping.")
                rmsd_list.append(-10)
                continue
            rmsd_list.append(rmsd)
        ranked_allowed_dihed_angles.append({
            "angle_indices": dihed_angle,
            "angle_atom_names": [transform_idx_to_name(angle_idx, mol) for angle_idx in dihed_angle], 
            "avg_rmsd": sum(rmsd_list)/len(rmsd_list)
            })
    ranked_allowed_dihed_angles = list(sorted(ranked_allowed_dihed_angles, key=lambda x: x["avg_rmsd"], reverse=True))
    non_duplicate_allowed_dihed_angles = []
    for dihed_angle in ranked_allowed_dihed_angles:
        should_append = True
        for non_duplicate in non_duplicate_allowed_dihed_angles:
            if abs(dihed_angle["avg_rmsd"] - non_duplicate["avg_rmsd"]) < 0.02:
                #TODO: better info messaging needed here, maybe include actual atom names? or skip altogether? not sure if this is useful
                log.info("Duplicate dihedral angle found when ranking dihedrals, skipping this dihedral angle.")
                should_append = False
                break
        if should_append:
            non_duplicate_allowed_dihed_angles.append(dihed_angle)
    num_unique_dihed_angles = len(non_duplicate_allowed_dihed_angles)
    log.info(f"Finished ranking all modifiable dihedral angles. {num_unique_dihed_angles} unique dihedral angles were identified and will be optimized in subsequent HAMR process.")
    return non_duplicate_allowed_dihed_angles



def transform_idx_to_name(atom_idx, mol):
    return mol.GetAtomWithIdx(atom_idx).GetMonomerInfo().GetName().strip()
    
def should_append_to_list(input_pdb, input_structs):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    ref_mol = Chem.MolFromPDBFile(input_pdb)
    for input_struct in input_structs:
        mol = Chem.MolFromPDBFile(input_struct["file_path"])
        try:
            if rdMolAlign.GetBestRMS(mol, ref_mol) < 0.005:
                return False
        except:
            continue
    return True
# def check_for_overlap(input_mol, threshold = 0):
#     from rdkit.Chem import rdGeometry
#     conf = input_mol.GetConformer()
#     for atom in input_mol.GetAtoms():
#         disallowed_idxs = [atom.GetIdx()]
#         atom_coords = conf.GetAtomPosition(atom.GetIdx())
#         for bond in atom.GetBonds():
#             disallowed_idxs.append(bond.GetOtherAtomIdx(atom.GetIdx()))
#         for inner_atom in input_mol.GetAtoms():
#             if inner_atom.GetIdx() in disallowed_idxs:
#                 continue
#             inner_atom_coords = conf.GetAtomPosition(inner_atom.GetIdx())
#             distance = rdGeometry.Point3D.Distance(atom_coords, inner_atom_coords)
#             if distance < threshold:
#                 return True
#     return False

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    print(find_amides(input_pdb="/Users/adam/Downloads/outputs_from_molec_replac/FINAL/PAR_ALPHA/PAR_NAT/TRIAL_0/ROUND_0/paritaprevir_alpha_0.pdb", template_smiles_string="Cc1cnc(cn1)C(=O)N[C@H]2CCCCC/C=C\\[C@@H]3C[C@]3(NC(=O)[C@@H]4C[C@H](CN4C2=O)Oc5c6ccccc6c7ccccc7n5)C(=O)NS(=O)(=O)C8CC8", log=log))
