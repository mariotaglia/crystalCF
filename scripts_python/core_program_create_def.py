import os
import shutil
import csv
import numpy as np
import subprocess
from collections import defaultdict
from transform_refs import calculate_center_ref, process_positions
from references.dependecies_init import list_reflexion
from function import run_command, read_DEF, write_DEF, path_carpeta, extract_params_init
from function import extract_R_bin, extract_references, update_particle_sizes, update_cdiva, extract_definitions
from function_part import generate_references_part_csv, extract_R_part, process_principal_part, process_secundario_part
from function_part import process_terciario_part
from create_defs import process_principal_binario, process_secundario_binario, process_terciario_binario, definitions_ref_edit

################### INICIO ##################
dir_inicial = os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python")

params_init = extract_params_init('init_params.txt')
name_bin = params_init['name']
n1 = params_init['n1']; n2 = params_init['n2']
R1_np = params_init['R1']

gamma_list = params_init['gamma list']
gamm_delta_dim = params_init['list gamma delta sum dim']
k_bin = params_init['num cell bin']
flag_reflexion = params_init["flag reflexion"]

delta_part = params_init["list delta part"]
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

# Not implemented yet
if flag_reflexion == True:
	n_k_list = list_reflexion(name_bin)
	n_k_bin = {"part1": n_k_list[0], "part2": n_k_list[1]}
else:
	n_k_bin = {"part1": n1*k_bin, "part2": n2*k_bin}
#
update_cdiva("DEFINITIONS.txt", name_bin)

for gamma_folder in gamma_folder_list:
	os.makedirs(f'gamma_{gamma_folder}', exist_ok=True)

for gamma_folder in gamma_folder_list:
	os.chdir(dir_inicial)
	DEF = os.path.join(dir_inicial, "DEFINITIONS.txt")
	tosubmit = os.path.join(dir_inicial, "tosubmit.sh")
	os.chdir(os.path.join(dir_inicial,f"gamma_{gamma_folder}"))

	dir_fuente = os.getcwd()
	folders = []
	folders = ["binary", "binary_ref/part1", "binary_ref/part2"]
	for part in ["part1", "part2"]:
		for cell in cell_part:
			folders.append(f"{part}/{cell}")
			folders.append(f"{part}/{cell}_ref")

	for folder in ["binary", "binary_ref", "part1", "part2"]:
		if os.path.exists(folder):
			try:
				#shutil.rmtree(folder)
				os.remove(os.path.join(folder, 'tot_references.csv'))
			except Exception as e:
				continue

	for folder in folders:
	    os.makedirs(folder, exist_ok=True)

	lines = read_DEF(DEF)
	gamma = float(gamma_folder.replace('_','.'))
	lines = update_particle_sizes(lines, gamma, R1_np, n_k_bin["part1"], n_k_bin["part2"])
	output_DEF = os.path.join("binary", "DEFINITIONS.txt")
	write_DEF(output_DEF, lines)

	os.chdir("binary")
	DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
	R1_np, R2_np = extract_R_bin(DEF, n1*k_bin)

	aL = float(run_command(f"python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {R2_np}"))
	delta_dim_bin = [entry for entry in gamm_delta_dim if entry["gamma"] == gamma]
	process_principal_binario(DEF, name_bin, delta_dim_bin, aL, n_k_bin, tosubmit, dir_fuente)
	os.chdir(dir_fuente)

	os.chdir("binary_ref")
	references_delta = extract_references("tot_references.csv")
	process_terciario_binario(os.getcwd(), name_bin, references_delta[1:], tosubmit, dir_fuente, n_k_bin)
	os.chdir(dir_fuente)
	print(f"created binary gamma = {gamma}")
	DEF_part = {}
	for cell in cell_part:
		DEF_part[cell] = os.path.join(dir_script,"references", f"DEFINITIONS_{cell}.txt")
	R_part = {"part1": R1_np, "part2": R2_np}
	
	for label in ["part1", "part2"]:
		for label_struc in cell_part:
			os.chdir(os.path.join(label, label_struc))
			dir_fuente_part = os.path.join(dir_fuente, label)
			shutil.copy(DEF_part[label_struc], os.path.join(os.getcwd(),"DEFINITIONS.txt"))
			DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
			lines = read_DEF(DEF)
			size_index = None
			for i, line in enumerate(lines):
				if line.strip() == "!particle semiaxis x y z in nm":
					size_index = i + 1
			for n in np.arange(0,k_part[label_struc]):
				lines[size_index + n] = f"{R_part[label]} {R_part[label]} {R_part[label]}\n"
			write_DEF(DEF, lines) 

			DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
			R_np = extract_R_part(DEF)
			aL = float(run_command(f"python3 {dir_script}/references/aL_min_{label_struc}.py {R_np}"))
			process_principal_part(DEF, delta_part[label_struc], aL, tosubmit, dir_fuente)
			os.chdir(dir_fuente)

	for label in ["part1", "part2"]:
	    for label_struc in cell_part:
	        DEF = os.path.join(dir_fuente, label, label_struc, "DEFINITIONS.txt")
	        os.chdir(os.path.join(dir_fuente, label, f"{label_struc}_ref"))
	        references_delta = extract_references("tot_references.csv")
	        process_terciario_part(os.getcwd(), references_delta[1:], DEF, tosubmit)
	        os.chdir(dir_fuente)

	print(f"created parts gamma = {gamma}")