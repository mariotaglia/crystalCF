import os
import subprocess
from function import extract_params_init

def editar_tosubmit(base_path, name_bin, gamma):
    tosubmit_path = os.path.join(base_path, 'tosubmit.sh')
    if os.path.exists(tosubmit_path):
        with open(tosubmit_path, 'r') as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("#SBATCH --job-name="):
                lines[i] = f"#SBATCH --job-name=\"{name_bin}_g_{gamma}\"\n"
                break
        with open(tosubmit_path, 'w') as file:
            file.writelines(lines)

def run_process_final(gamma_list, name_bin):
    dir_origin = os.getcwd()
    for i, gamma in enumerate(gamma_list):
        base_path = os.path.join(dir_origin,f"gamma_{gamma}")
        paths = []
        
        if i > 0:
            search_path = [d for d in os.walk(base_path) if f"gamma_{gamma}/part1" not in d[0]]
        else:
            search_path = os.walk(base_path)

        for root, _, files in directorios_a_buscar:
            if 'tosubmit.sh' in files:
                print(root)
                paths.append(root)

        for dir1 in paths:
            editar_tosubmit(dir1, name_bin, gamma)
            os.chdir(dir1)
            os.system("sbatch tosubmit.sh")

################### START ##################
dir_origin = os.getcwd()
dir_script = os.path.join(dir_origin,'modulos')

params_init = extract_params_init('init_params.txt')
name_bin = params_init['name']
gamma_list = params_init['gamma list']
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

os.chdir(dir_origin)
run_process_final(gamma_folder_list, name_bin)