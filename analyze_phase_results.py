import hamr_calc_struct_fact

def main(input_dir, input_mtz, intensity_column, sigma_column):
    import os
    import logging
    from pathlib import Path
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    solutions = []

    for file in os.listdir(input_dir):
        #TODO: is this appropriate? how to use customizability for RMSD calculation debugging?
        if ".sol" not in file:
                continue
        file_path = f"{input_dir}/{file}"
        with open(file_path, "r") as f:
                solutions.append(extract_solu_stats(
                    solu_string=f.read(), 
                    model_name=file.split(".sol")[0], 
                    log=log,
                    file_path=file_path, 
                    input_mtz=input_mtz,
                    intensity_column=intensity_column,
                    sigma_column=sigma_column
                    ))
    solutions = list(sorted(solutions, key = lambda x : x["llg"], reverse=True))
    return solutions

def extract_solu_stats(solu_string, model_name, log, file_path, input_mtz, intensity_column, sigma_column):
    import re
    fixed_file_path = re.sub(".sol",".1.pdb", file_path)
    extracted_data = {"name": model_name, "file_path": fixed_file_path, "llg": 0.0, "tfz": 0.0}
    solu_string = solu_string.strip()
    try:
        tfz_scores = re.findall(r"TFZ==[0-9]*\.[0-9]*", solu_string)
        extracted_data["tfz"] = float(tfz_scores[-1].split("==")[-1])
    except:
        log.warn(f"Failed extracting TFZ value from .sol PHASER output for {model_name}. Continuing.")
    try:
        llg_scores = re.findall(r"LLG=[0-9]*\.?[0-9]*", solu_string)
        extracted_data["llg"] = float(llg_scores[-1].split("=")[-1])
    except:
        log.warn(f"Failed extrating LLG value from .sol PHASER output for {model_name}. Continuing.")
    try:
        space_groups = re.findall(r"SOLU SPAC .*", solu_string)
        extracted_data["space_group"] = re.sub(r"SOLU SPAC|SOLU.*", "", space_groups[-1]).strip()
    except:
        log.warn(f"Failed extracting space group from .sol PHASER output for {model_name}. Continuing.")
    try:
         r_factor = hamr_calc_struct_fact.main(input_pdb=fixed_file_path, input_mtz=input_mtz, intensity_column=intensity_column, sigma_column=sigma_column)
         extracted_data["r_factor"] = r_factor
    except:
         log.warn(f"Failed calculating R-factor for PHASER output for {file_path}. Setting to 1 and continuing.")
         extracted_data["r_factor"] = 1
    return extracted_data

if __name__ == "__main__":
    results = main("/Users/adam/Downloads/outputs_from_molec_replac/FINAL/PAR_ALPHA/PAR_12A/CONFORGE/PHASER", input_mtz="/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.mtz", intensity_column="I", sigma_column="SIGI")
    llgs = [x["llg"] for x in results]
    tfzs = [x["tfz"] for x in results]
    # for result in results:
    #      result["r_fac"] = hamr_calc_struct_fact.main("/Users/adam/Downloads/inputs_for_molec_replac/paritaprevir_alpha.mtzz", result["file_path"], "I", "SIGI")
    for entry in list(sorted(results, key=lambda x: x["r_factor"], reverse=False))[:20]:
         print(f"{entry['name']} -- llg: {entry['llg']}, r_factor: {entry['r_factor']}")
    print(sum(llgs)/len(llgs))
    print(sum(tfzs)/len(tfzs))
        