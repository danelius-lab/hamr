import logging

def main(input_dir):
    import os
    log = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    log.info("Extracting results after PHENIX refinement.")
    results = []
    for f in os.listdir(input_dir):
        if ".log" in f:
            r_facs = extract_r_factor(f"{input_dir}/{f}")
            results.append({"r_work": r_facs[0], "r_free": r_facs[1], "name": f})
    results = list(sorted(results, key = lambda x: x["r_free"]))
    log.info("Results from this complete HAMR trial after refinement are below:")
    for result in results:
        log.info(f"{result['name']} -- R_free: {result['r_free']}, R_work: {result['r_work']}")
    return results[0]["r_free"] < 0.3 and results[0]["r_work"] < 0.3 and results[0]["r_free"] - results[0]["r_work"] < 0.05



def extract_r_factor(input_path):
    import re
    import os
    for f in os.listdir(input_path):
        if ".log" not in f: 
            continue    
        with open(f"{input_path}/{f}", "r") as inner_f:
            log_string = inner_f.read()
            r_factor_regex = re.compile(r"Final R-work = [0-9].[0-9]*, R-free = [0-9].[0-9]*")
            res = re.findall(r_factor_regex, log_string)
            if len(res) < 1:
                return [1, 1]
            split_res = res[0].split(",")
            r_factors = []
            for split_str in split_res:
                r_factors.append(float(split_str.split("=")[-1]))
            inner_f.close()
            return r_factors

