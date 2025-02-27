import sys
import json
import logging
from pathlib import Path

import conforge_create_confs
import phase
import analyze_phase_results
def main(settings, log):
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

    output_conformer_path = conforge_create_confs.main(
        input_smiles_path=input_smiles_path, 
        output_path=output_path, 
        conformer_prefix=conformer_prefix, 
        max_time=max_time,
        max_confs=max_confs,
        e_window=e_window,
        min_rmsd=min_rmsd
        )
    try:
        input_mtz = settings["input_mtz"]
        output_path = settings["output_path"]
    except:
        log.exception("Necessary settings (input_mtz, output_pdb) are missing or invalid in JSON settings. Please provide theses field in the correct format. Exiting.")
        exit(1)
    try:
        should_log_phaser = settings["should_log_phaser"]
    except:
        log.warning("Failed extrating optional setting (should_log_phaser) from settings JSON. Using default setting.")
        should_log_phaser = False
    phase.main(
        input_dir=output_conformer_path,
        should_log_phaser = should_log_phaser,
        input_mtz=input_mtz,
        output_path=output_path
    )
    try:
        solutions = analyze_phase_results.main(output_path)
        log.info("Initial MR solutions from CONFORGE generated conformers")
        print_phaser_results(phaser_results = solutions, log=log)
    except:
        log.exception("Failed extracting solutions from PHASER output. Exiting.")
        exit(1)
    current_soln = solutions[0]["file_path"]

    return


def print_phaser_results(phaser_results, log):
    for entry in phaser_results:
        log.info(f"name: {entry["name"]}, tfz: {entry["tfz"]}, llg: {entry["llg"]}, space group: {entry["space_group"]}")

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    # TODO: add in more information for startup, email, paper, etc...
    log.info("Starting High-throughput Automated Replacement (HAMR)")
    if len(sys.argv) < 2:
        log.exception("No settings provided by user. Exiting.")
        exit(1)
    settings_path = sys.argv[2]
    with open(settings_path, "r") as f:
        try:
            settings = json.load(f)
            f.close()
        except:
            log.exception("Invalid JSON settings, either incorrect formatting or invalid file path provided. Exiting.")
            exit(1)
    main(settings=settings, log=log)


        
    

