import os
import subprocess
import time

def run_process_final(gamma_list, name_bin, run_program_file):
    dir_origin = os.getcwd()
    
    base_path = os.path.join(dir_origin,f"sim_part1")
    paths = []
        
    search_path = os.walk(base_path)

    for root, _, files in search_path:
        if 'fitpairL12.dat' in files:
            paths.append(root)

    for dir1 in paths:
        os.chdir(dir1)
        os.system(f"python3 {run_program_file}")
        os.system("wait")

    for i, gamma in enumerate(gamma_list):
        base_path = os.path.join(dir_origin,f"gamma_{gamma}","binary")
        paths = []
        
        search_path = os.walk(base_path)
        for root, _, files in search_path:
            if 'fitpairL12.dat' in files:
                #paths.append(root)
                print(root)

        for dir1 in paths:
            os.chdir(dir1)
            os.system(f"python3 {run_program_file}")
            os.system("wait")

    for i, gamma in enumerate(gamma_list):
        base_path = os.path.join(dir_origin,f"gamma_{gamma}","part2")
        paths = []
        
        search_path = os.walk(base_path)
        for root, _, files in search_path:
            if 'fitpairL12.dat' in files:
                #continue
                paths.append(root)

        for dir1 in paths:
            os.chdir(dir1)
            os.system(f"python3 {run_program_file}")
            os.system("wait")

def read_DEF(file_path):
    """Extract the lines from DEF."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines

def extract_params_init(params_init):
    data = {
        "name": None,
        "gamma list": []
    }

    lines = read_DEF(params_init)
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if line == "!name":
            data["name"] = lines[i + 1].strip("\n")
            i += 1
        elif line == "!list gamma":
            data["gamma list"] = [float(x) for x in lines[i+1].strip("[]\n").split(",")]
            i += 1
        i += 1
    return data

################### START ##################
dir_origin = os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python/pairwise_additive_F")
run_program_file = os.path.join(dir_script,"calc_F_pair.py")

params_init = extract_params_init('init_params.txt')
name_bin = params_init['name']
gamma_list = params_init['gamma list']
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]
os.chdir(dir_origin)
run_process_final(gamma_folder_list, name_bin, run_program_file)