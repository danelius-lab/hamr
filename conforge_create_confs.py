from sys import exit

def generateConformers(mol, conf_gen):
     from CDPL import ConfGen
     ConfGen.prepareForConformerGeneration(mol)
     status = conf_gen.generate(mol)
     num_confs = conf_gen.getNumConformers()
     if status == ConfGen.ReturnCode.SUCCESS or status == ConfGen.ReturnCode.TOO_MUCH_SYMMETRY:
         conf_gen.setConformers(mol)
     else:
         num_confs = 0
     return (status, num_confs)

def split_sdf_to_pdb(log, sdf_input_path, pdb_output_path, conformer_prefix):
    from rdkit import Chem
    try:
        mols = Chem.SDMolSupplier(sdf_input_path, removeHs=False)
        count = 0
        for mol in mols:
            if mol == None:
                log.warning(f"Found invalid molecule in output CONFORGE conformer SDF at position {count}. Continuing.")
                continue
            Chem.MolToPDBFile(mol, f"{pdb_output_path}/{conformer_prefix}_{count}.pdb")
            count += 1
        return pdb_output_path
    except:
        log.exception(f"Failed splitting generated conformers from SDF to PDB format. Exiting.")
        exit(1)

def main(input_smiles_path, output_path, conformer_prefix, max_time=7200, min_rmsd=0.5, e_window=20, max_confs=100):
    from CDPL import Chem, ConfGen
    import logging
    import os
    
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    log.info("Starting initial conformer generation with CONFORGE. This may take several minutes without updating.")
    os.makedirs(output_path, exist_ok=True)
    # return output_path
    # output_sdf_path = f"{output_path}/{conformer_prefix}.sdf"
    # return split_sdf_to_pdb(sdf_input_path=output_sdf_path, pdb_output_path=output_path, conformer_prefix=conformer_prefix, log=log)
    try:
        reader = Chem.MoleculeReader(input_smiles_path)
    except:
        log.exception(f"Invalid or missing smiles file for initial conformer generation at: {input_smiles_path}")
        exit(1)
    mol = Chem.BasicMolecule()
    # TODO: update this with more settings (namely number of sampled conformers disabled limiting, sampling mode)
    try:
        conf_gen = ConfGen.ConformerGenerator()
        conf_gen.settings.timeout = max_time * 1000
        conf_gen.settings.minRMSD = min_rmsd
        conf_gen.settings.energyWindow = e_window
        conf_gen.settings.maxNumOutputConformers = max_confs
        conf_gen.settings.genCoordsFromScratch = True
    except:
        log.exception(f"Error initializing CONFORGE conformer generator with the following settings: \n time_limit (s):{max_time}\n minimum_rmsd (Å): {min_rmsd}\n energy_window (kcal/mol): {e_window}\n max_num_conformers: {max_confs}")
        exit(1)
    reader.read(mol)
    try:
        _, num_confs = generateConformers(mol=mol, conf_gen=conf_gen)
    except Exception as e:
        # TODO: add better error reporting here
        log.exception(f"Failed initial CONFORGE conformer generation. Exiting. Detailed exception information:\n {e}")
        exit(1)
    log.info(f"Succesfully generated {num_confs} conformers with CONFORGE.")
    try:
        output_sdf_path = f"{output_path}/{conformer_prefix}.sdf"
        writer = Chem.MolecularGraphWriter(output_sdf_path)
        writer.write(mol)
    except Exception as e:
        log.exception(f"Failed writing conformers to specified output path: {output_path}. Exiting. Detailed exception information:\n {e}")
        exit(1)
    return split_sdf_to_pdb(sdf_input_path=output_sdf_path, pdb_output_path=output_path, conformer_prefix=conformer_prefix, log=log)

    
    
