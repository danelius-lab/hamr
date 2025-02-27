import os
import time
import subprocess
import logging
def main(input_dir, input_mtz, output_path, should_log_phaser=False):
    log = logging.getLogger(__name__)
    for file in list(sorted(os.listdir(input_dir))):
        if ".pdb" not in file:
            continue
        #TODO: add identity logic for MR
        prefix = file.split("/")[-1].split(".pdb")[0]
        run_phaser(
            input_mtz=input_mtz,
            input_pdb=file,
            output_path=output_path,
            prefix=prefix,
            log=log,
            should_log_phaser=should_log_phaser
            )


def run_phaser(input_mtz, input_pdb, output_path, prefix, log, should_log_phaser=False, identity=1):
    phaser_exec = os.environ.get("HAMR_PHASER_EXEC")
    if not phaser_exec:
        log.exception("Could not retrieve PHASER executable from environment variables. Please specify the location of your PHASER executable at $HAMR_PHASER_EXEC (e.g. export HAMR_PHASER_EXEC=/path/to/your/executable). Exiting.")
        exit(1)
    #COMPOSITION ATOM H NUMBER 1 is a workaround for compositions of small molecules not fitting in unit cell volume
    #the composition of the ASU is automatically correct by PHASER without failure (most of the time)
    #TODO: make this more customizable (COMPOSITION, SEARCH METHOD, PURGE ROT/TRA/RNP lists)
    phaser_str = f""" << eof
    MODE MR_AUTO
    HKLIN {input_mtz}
    COMPOSITION ATOM H NUMBER 1
    ENSEMBLE HAMR_SEARCH PDBFILE {input_pdb} IDENTITY {identity}
    ENSEMBLE HAMR_SEARCH DISABLE CHECK ON
    ENSEMBLE HAMR_SEARCH HETATM ON
    FORMFACTORS ELECTRON
    PACK SELECT ALL
    ELLG TARGET 225
    SGALTERNATIVE SELECT ALL
    XYZOUT ON ENSEMBLE ON
    XYZOUT ON PACKING ON 
    TOPFILES 1
    KEYWORDS ON
    ZSCORE USE OFF
    ROOT {output_path}/{prefix}_PHASER
    SEARCH ENSEMBLE HAMR_SEARCH
    SEARCH METHOD FAST
    eof"""
    try:
        log.info(f"Starting PHASER for {prefix}")
        start = time.time()
        proc = subprocess.run([phaser_exec + phaser_str], shell=True, capture_output=True, close_fds=False, start_new_session=True)
        end = time.time()
        total_time = round(end - start, 1)
        log.info(f"Finished PHASER for {prefix} in {total_time} seconds")
        if should_log_phaser:
            phaser_output = proc.stdout.decode()
            with open(f"{output_path}/phaser_stdout.txt", "w") as f:
                f.write(phaser_output)
                f.close()
            log.info(phaser_output)
    except:
        log.warn(f"Failed PHASER for {prefix}. Continuing with PHASER runs.")
    
    

