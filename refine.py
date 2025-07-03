import logging
from sys import exit

import hamr_create_confs
import analyze_refine_results
def main(input_pdbs, 
        input_mtz, 
        restraint_cif,
        num_refine_cycles, 
        r_free_fraction,
        refinement_columns,
        output_dir):
    import os
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    log.info(f"Starting PHENIX refinement process for top {len(input_pdbs)} HAMR solutions.")
    count = 0
    r_facs = []
    for input_pdb in input_pdbs:
        
        pdb_header = hamr_create_confs.extract_pdb_header(input_pdb["file_path"], log)
        space_group_string = pdb_header["space_group"]
        unit_cell_dimensions = pdb_header["unit_cell_dimensions"]
        prefix = input_pdb["file_path"].split("/")[-1].split(".1.pdb")[0]
        output_refine_dir = f"{output_dir}/{prefix}"
        os.makedirs(output_refine_dir, exist_ok=True)
        r_fac = refine(
            input_pdb=input_pdb["file_path"],
            input_mtz=input_mtz,
            restraint_cif=restraint_cif,
            space_group_string=space_group_string,
            unit_cell_dimensions=unit_cell_dimensions,
            r_free_fraction=r_free_fraction,
            refinement_columns=refinement_columns,
            output_dir=output_refine_dir,
            num_refine_cycles=num_refine_cycles,
            count=count,
            log=log
        )
        count += 1
        r_facs.append(r_fac)
    results = list(sorted(r_facs, key = lambda x: x["r_free"]))
    log.info("Results from this complete HAMR trial after refinement are below:")
    str_result_list = []
    for result in results:
        str_res = f"{result['name']} -- R_work: {result['r_work']}, R_free: {result['r_free']}"
        log.info(str_res)
        str_result_list.append(str_res)
    with open(f"{output_dir}/refine_summary.txt", "w") as f:
        f.write("\n".join(str_result_list))
        f.close()
    return results[0]["r_free"] < 0.3 and results[0]["r_work"] < 0.3 and results[0]["r_free"] - results[0]["r_work"] < 0.05

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
    r_factors = analyze_refine_results.extract_r_factor(output_dir)
    r_free = r_factors[1]
    r_work = r_factors[0]
    log.info(f'Finished PHENIX refinement for {input_pdb} -- R_work: {r_work}, R_free: {r_free}')
    del proc
    return {"r_free": r_free, "r_work": r_work, "name": input_pdb}

