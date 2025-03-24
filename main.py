import sys
import json
import logging
import os
from pathlib import Path
from sys import exit

import conforge_create_confs
import phase
import analyze_phase_results
import hamr
def main(settings, log):
    phaser_exec = os.environ.get("HAMR_PHASER_EXEC")
    phenix_setup_script = os.environ.get("HAMR_PHENIX_SETUP_SH")
    if phaser_exec == None:
        log.exception("Could not retrieve PHASER executable from environment variables. Please specify the location of your PHASER executable at $HAMR_PHASER_EXEC (e.g. export HAMR_PHASER_EXEC=/path/to/your/executable). Exiting.")
        exit(1)
    if phenix_setup_script == None:
        log.exception("Could not retrieve PHENIX setup script from environment variables. Please specify the location of your PHENIX setup script at $HAMR_PHENIX_SETUP_SH (e.g. export HAMR_PHENIX_SETUP_SH=/path/to/your/phenix_setup.sh). Exiting.")
        exit(1)
    try:
        input_smiles_path = settings["input_smiles_path"]
        output_path = settings["output_path"]
        conformer_prefix = settings["conformer_prefix"]
    except:
        log.exception("Necessary settings (input_smiles_path, output_path, and/or conformer_prefix) missing or invalid in JSON settings. Please provide these fields in the correct format. Exiting.")
        exit(1)
    try:
        Path(output_path).mkdir(exist_ok=True, parents=True)
    except:
        log.exception("Failed creating or viewing the output directory specified in JSON settings. Ensure that this is a proper path. Exiting.")
        exit(1)
    try:
        max_time = settings["conforge_max_time"]
        min_rmsd = settings["conforge_min_rmsd"]
        e_window = settings["conforge_energy_window"]
        max_confs = settings["conforge_max_num_confs"]
    except:
        log.info("Failed extracting optional settings (conforge_max_time, conforge_min_rmsd, conforge_energy_window, conforge_max_num_confs) from settings JSON. Using default settings.")
        max_time = 7200
        max_confs = 100
        e_window = 20
        min_rmsd = 0.5
    os.makedirs(output_path, exist_ok=True)
    output_conformer_path = conforge_create_confs.main(
        input_smiles_path=input_smiles_path, 
        output_path=f"{output_path}/CONFORGE", 
        conformer_prefix=conformer_prefix, 
        max_time=max_time,
        max_confs=max_confs,
        e_window=e_window,
        min_rmsd=min_rmsd
        )
    try:
        input_mtz = settings["input_mtz"]
        output_path = settings["output_path"]
        r_factor_columns = settings["r_factor_columns"]
    except:
        log.exception("Necessary settings (input_mtz, output_path, r_factor_columns) are missing or invalid in JSON settings. Please provide theses field in the correct format. Exiting.")
        exit(1)
    try:
        should_log_phaser = settings["should_log_phaser"]
    except:
        log.warning("Failed extrating optional setting (should_log_phaser) from settings JSON. Using default setting.")
        should_log_phaser = False
    phase.main(
        input_dir=output_conformer_path,
        should_log_phaser=should_log_phaser,
        input_mtz=input_mtz,
        output_path=f"{output_conformer_path}/PHASER"
    )
    try:
        solutions = analyze_phase_results.main(
            input_dir=f"{output_conformer_path}/PHASER",
            input_mtz=input_mtz,
            intensity_column=r_factor_columns[0],
            sigma_column=r_factor_columns[1])
        log.info("Initial MR solutions from CONFORGE generated conformers")
        print_phaser_results(phaser_results=solutions, log=log)
    except:
        log.exception("Failed extracting solutions from PHASER output. Exiting.")
        exit(1)
    #TODO: consider error handling? idk if this is necessary here
    
    current_soln = solutions[0]["file_path"]
    current_trial = 0
    has_been_solved = False
    while not has_been_solved and current_trial < len(solutions):
        # if "268" not in current_soln:
        #     current_trial += 1
        #     current_soln = solutions[current_trial]["file_path"]
        #     continue
        try:
            restraint_cif = settings["restraint_cif"]
            smiles_string = settings["smiles_string"]
            refinement_columns = settings["refinement_columns"]            

        except:
            log.exception("Necessary settings (restraint_cif, smiles_string, refinement_column) are missing or invalid in JSON settings. Please provide theses field in the correct format. Exiting.")
            exit(1)
        try:
            num_refine_cycles = settings["num_refine_cycles"]
            angle_step = settings["angle_step"]
            num_to_persist = settings["num_to_persist"]
            faulty_conformer_rmsd_cutoff = settings["faulty_conformer_rmsd_cutoff"]
            r_free_fraction = settings["r_free_fraction"]
            should_force_trans_amides = settings["should_force_trans_amides"]
            fraction_to_phase = settings["fraction_to_phase"]
        except:
            log.warning("Failed extracting optional setting (num_refine_cycles, angle_stop, num_to_persist, faulty_conformer_rmsd_cutoff, r_free_fraction, should_force_trans_amides) from settings JSON. Using default setting.")
            num_refine_cycles = "5"
            angle_step = 10
            num_to_persist = 15
            faulty_conformer_rmsd_cutoff = 0
            r_free_fraction = "0.05"    
            should_force_trans_amides = True
            fraction_to_phase = 0.25
        
        has_been_solved = hamr.main(
            input_pdb=current_soln,
            input_mtz=input_mtz,
            restraint_cif=restraint_cif,
            output_path=f"{output_path}/TRIAL_{current_trial}",
            smiles_string=smiles_string,
            r_free_fraction=r_free_fraction,
            refinement_columns=refinement_columns,
            r_factor_columns=r_factor_columns,
            num_refine_cycles=num_refine_cycles,
            num_to_persist=num_to_persist,
            prefix=conformer_prefix,
            angle_step=angle_step,
            faulty_conformer_rmsd_cutoff=faulty_conformer_rmsd_cutoff,
            should_force_trans_amides=should_force_trans_amides,
            fraction_to_phase=fraction_to_phase
        )
        current_trial += 1
        if not has_been_solved and current_trial < len(solutions):
            current_soln = solutions[current_trial]["file_path"]
    if current_trial < len(solutions):
        log.info(f"Succesfully found a solution, the data and results for this is stored in: {output_path}/TRIAL_{current_trial}. Exiting.")
        exit(0)
    else:
        log.info(f"Failed to identify a solution after {len(solutions)} trials. Bummer it didn't work out. Exiting.")
        exit(1)

def print_phaser_results(phaser_results, log):
    for entry in phaser_results:
        name = entry["name"]
        tfz = entry["tfz"]
        llg = entry["llg"]
        space_group = entry["space_group"]
        log.info(f"name: {name}, tfz: {tfz}, llg: {llg}, space group: {space_group}")

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    # TODO: add in more information for startup, email, paper, etc...
    log.info("Starting High-throughput Automated Replacement (HAMR)")
    if len(sys.argv) < 2:
        log.exception("No settings provided by user. Exiting.")
        exit(1)
    settings_path = sys.argv[-1]
    with open(settings_path, "r") as f:
        try:
            settings = json.load(f)
            f.close()
        except:
            log.exception("Invalid JSON settings, either incorrect formatting or invalid file path provided. Exiting.")
            exit(1)
    main(settings=settings, log=log)


        
    

