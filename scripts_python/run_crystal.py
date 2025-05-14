import os
import subprocess
import time

def editar_tosubmit(base_path, name_bin, nseg, cov):
    tosubmit_path = os.path.join(base_path, 'tosubmit.sh')
    if os.path.exists(tosubmit_path):
        with open(tosubmit_path, 'r') as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("#SBATCH --job-name="):
                lines[i] = f"#SBATCH --job-name=\"{name_bin}_n_{nseg}_c_{cov}\"\n"
                break
        with open(tosubmit_path, 'w') as file:
            file.writelines(lines)

def run_process_final(nseg_folder_list, cov_folder_list, name_bin):
    dir_origin = os.getcwd()

    for nseg in nseg_folder_list:
        for cov in cov_folder_list:
            base_path = os.path.join(dir_origin,f"nseg_{nseg}",f"cov_{cov}")
            paths = []
            
            search_path = os.walk(base_path)

            for root, _, files in search_path:
                if 'tosubmit.sh' in files:
                    paths.append(root)

            for dir1 in paths:
                editar_tosubmit(dir1, name_bin, nseg, cov)
                os.chdir(dir1)
                os.system("sbatch tosubmit.sh")
                time.sleep(0.01)

def read_DEF(file_path):
    """Extract the lines from DEF."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines

def extract_params_init(params_init):
    data = {
        "name": None,
        "nseg list": [],
        "coverage list": []
    }

    lines = read_DEF(params_init)
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if line == "!name":
            data["name"] = lines[i + 1].strip("\n")
            i += 1
        elif line == "!list long":
            data["nseg list"] = [float(x) for x in lines[i+1].strip("[]\n").split(",")]
            i += 1
        elif line == "!list coverage":
            data["coverage list"] = [float(x) for x in lines[i+1].strip("[]\n").split(",")]
            i += 1
        i += 1
    return data

################### START ##################
dir_origin = os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python")

params_init = extract_params_init('init_params.txt')
name_bin = params_init['name']
Nseg_list = params_init['nseg list']
Cov_list = params_init['coverage list']
nseg_folder_list = ["{:.0f}".format(n).replace('.','_') for n in Nseg_list]
cov_folder_list = ["{:.2f}".format(c).replace('.','_') for c in Cov_list]

os.chdir(dir_origin)
run_process_final(nseg_folder_list, cov_folder_list, name_bin)