import logging

import hamr_create_confs
def main(input_pdbs, 
        input_mtz, 
        restraint_cif,
        num_refine_cycles, 
        r_free_fraction,
        refinement_columns,
        output_dir):
    log = logging.GetLogger(__name__)
    log.info(f"Starting PHENIX refinement process for top {len(input_pdbs)} HAMR solutions.")
    count = 0
    for input_pdb in input_pdbs:
        count += 1
        pdb_header = hamr_create_confs.extract_pdb_header(input_pdb, log)
        space_group_string = pdb_header["space_group"]
        unit_cell_dimensions = pdb_header["unit_cell_dimensions"]
        refine(
            input_pdb=input_pdb,
            input_mtz=input_mtz,
            restraint_cif=restraint_cif,
            space_group_string=space_group_string,
            unit_cell_dimensions=unit_cell_dimensions,
            r_free_fraction=r_free_fraction,
            refinement_columns=refinement_columns,
            output_dir=output_dir,
            num_refine_cycles=num_refine_cycles,
            count=count,
            log=log
        )

def shell_source(script, log):
    import subprocess
    import os
    
    log.info("Setting up PHENIX environment variables.")
    pipe = subprocess.Popen(f"/bin/bash {script}", stdout=subprocess.PIPE, shell=True)
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines()))
    os.environ.update(env)

def refine(
          input_pdb,
          input_mtz,
          restraint_cif,
          space_group_string,
          unit_cell_dimensions,
          r_free_fraction,
          refinement_columns,
          output_dir,
          num_refine_cycles,
          count,
          log
          ):
    import os
    import subprocess
    log.info(f"Starting PHENIX refinement for {input_pdb} ")
    phenix_setup_script = os.environ.get("HAMR_PHENIX_SETUP_SH")

    with open("./refinement_params","r") as f:
        param_str = f.read()
        f.close()
    param_str = param_str.replace("REPLACE_PDB", input_pdb)
    param_str = param_str.replace("REPLACE_MTZ", input_mtz)
    param_str = param_str.replace("REPLACE_RESTRAINT", restraint_cif)
    param_str = param_str.replace("REPLACE_SPACEGROUP", space_group_string)
    param_str = param_str.replace("REPLACE_UNIT_CELL", unit_cell_dimensions)
    param_str = param_str.replace("REPLACE_RFLAG_FRACTION", r_free_fraction)
    param_str = param_str.replace("REPLACE_COLUMN", refinement_columns)
    param_str = param_str.replace("REPLACE_NUM_CYCLES", num_refine_cycles)
    with open(f"{output_dir}/phenix_params_{count}","w") as f:
        f.write(param_str)
        f.close()
                
    shell_source(phenix_setup_script, log)
    proc = subprocess.run([f"cd {output_dir} && source {phenix_setup_script} && phenix.refine ./phenix_params_{count} overwrite=true refinement.input.symmetry_safety_check=warning" ], capture_output=True, shell=True)
    log.info(f'Finished PHENIX refinement for {input_pdb}')
    
    
                