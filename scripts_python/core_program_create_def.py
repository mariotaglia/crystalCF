import os
import shutil
import csv
import numpy as np
import subprocess
from collections import defaultdict
from transform_refs import calculate_center_ref, process_positions
from function import run_command, read_DEF, write_DEF, path_carpeta, extract_params_init
from function import extract_R_bin, extract_references, update_particle_sizes, extract_definitions
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
delta_bin = params_init['list delta bin']
dims_sum_bin = params_init['list sum dim bin']
k_bin = params_init['num cell bin']

delta_part = params_init["list delta part"]
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

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
	lines = update_particle_sizes(lines, gamma, R1_np, n1*k_bin, n2*k_bin)
	output_DEF = os.path.join("binary", "DEFINITIONS.txt")
	write_DEF(output_DEF, lines)
	print(f"Archivo {output_DEF} generado correctamente con gamma = {gamma}")

	os.chdir("binary")
	DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
	R1_np, R2_np = extract_R_bin(DEF, n1*k_bin)

	aL = float(run_command(f"python3 {dir_script}/references/aL_min_{name_bin}.py {R1_np} {R2_np}"))
	n = {"part1": n1, "part2": n2}
	process_principal_binario(DEF, delta_bin, aL, n, k_bin, tosubmit, dir_fuente, dims_sum_bin[gamma])
	os.chdir(dir_fuente)

	os.chdir("binary_ref")
	references_delta = extract_references("tot_references.csv")
	process_terciario_binario(os.getcwd(), references_delta[1:], tosubmit, dir_fuente, n, k_bin)
	os.chdir(dir_fuente)

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

