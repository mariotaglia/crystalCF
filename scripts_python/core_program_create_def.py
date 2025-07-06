import os
import sys
import shutil
import csv
import numpy as np
import pandas as pd
import subprocess
from collections import defaultdict
from transform_refs import calculate_center_ref, process_positions
from references.dependecies_init import list_reflexion
from function import run_command, read_DEF, write_DEF, path_carpeta, extract_params_init
from function import extract_R_bin, extract_references, update_particle_sizes, update_R1, extract_definitions, gamma_calc
from cdiva_function import update_cdiva
from function_part import generate_references_part_csv, extract_R_part, process_principal_part, process_secundario_part
from function_part import process_terciario_part
from create_defs import process_principal_binario, process_secundario_binario, process_terciario_binario, definitions_ref_edit

################### INICIO ##################
dir_inicial = os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python")
params_init = extract_params_init('init_params.txt', False)
name_bin = params_init['name']
n1 = params_init['n1']; n2 = params_init['n2']
R1_np = params_init['R1']

gamma_list = params_init['gamma list']
gamm_delta_dim = params_init['list gamma delta sum dim']
k_bin = params_init['num cell bin']
n_k_bin = {"part1": n1*k_bin, "part2": n2*k_bin}
flag_reflexion = params_init["flag reflexion binary"]
flag_reflexion_part = params_init["flag reflexion part"]

part_delta_dim = params_init['list part delta sum dim']
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

update_R1("DEFINITIONS.txt", n_k_bin["part1"], R1_np)
DEF = os.path.join(dir_inicial, "DEFINITIONS.txt")
lines = read_DEF(DEF)
for i, line in enumerate(lines):
	if line.strip() == "!seed":
		size_index = i + 1
		seed = lines[size_index]
		seed_lig = lines[size_index + 1]
	elif line.strip() == "!properties of ligand chains":
		size_index = i + 1
		nseg = lines[size_index].split()[1]

n1_k_bin = n_k_bin["part1"]; n2_k_bin = n_k_bin["part2"]
sections_info = [
    ("! Rotation", 3, n1_k_bin+n2_k_bin),
    ("! coverage", 1, n1_k_bin+n2_k_bin),
    ("!Surface-polymer atraction", 1, n1_k_bin+n2_k_bin),
    ("!chains lenght", 1, n1_k_bin+n2_k_bin)
]

sections_found = []
for key, lines_per_particle, tot_particles in sections_info:
    for i, line in enumerate(lines):
        if line.strip() == key:
            start_index = i + 1
            sections_found.append((key, start_index, lines_per_particle, tot_particles))
            break
    else:
        print(f"Advertencia: No se encontró la sección {key}.")

configs = [("part1", 0, n1_k_bin), ("part2", n1_k_bin, n2_k_bin)]  # (Name, Offset, number of particles)
cov = {"part1": None, "part2": None}
chain_lenght = {"part1": None, "part2": None}
for label, offset, num_particles in configs:
    modified_lines = lines.copy()

    for key, start_index, lines_per_particle, tot_particles in sorted(sections_found, key=lambda x: x[1], reverse=True):
        block_length = tot_particles * lines_per_particle
        start_offset = start_index + offset*lines_per_particle
        end_offset = start_offset + num_particles*lines_per_particle
        # Extract the particle params from DEFINITIONS
        new_block = modified_lines[start_offset:end_offset]
        if key == "! coverage":
        	cov[label] = [elem.replace('\n', '') for elem in new_block][0]
        elif key == "!chains lenght":
        	chain_lenght[label] = [elem.replace('\n', '') for elem in new_block][0]

for label in ['part1','part2']:
	if chain_lenght[label] == None:
		chain_lenght[label] = nseg

k_aL = {"kx": 1,"ky": 1,"kz": 1}
if flag_reflexion == True:
	PBC = []
	for i, line in enumerate(lines):
		if line.strip() == "!PBC PBC xmin xmax ymin ymax zmin zmax, 1=yes, 2=wall, 0=bulk":
			PBC.append([float(x) for x in lines[i+1].split()[1:]])
			break
	PBC = PBC[0]
	if PBC[0] != PBC[1] or PBC[2] != PBC[3] or PBC[4] != PBC[5]:
		print("PBC bad entry error.")
		sys.exit()
	else:
		for i, k in enumerate(k_aL):
			if PBC[2*i] == 3:
				k_aL[k] = 2
			else:
				k_aL[k] = 1

for gamma_folder in gamma_folder_list:
	os.makedirs(f'gamma_{gamma_folder}', exist_ok=True)

part1_folder = ["sim_part1/binary_ref","sim_part1/binary_ref/part1"]
for cell in cell_part:
	part1_folder.append(f"sim_part1/{cell}")
	part1_folder.append(f"sim_part1/{cell}_ref")

for folder in part1_folder:
	os.makedirs(folder, exist_ok=True)
	if os.path.exists(folder):
		try:
			#shutil.rmtree(folder)
			os.remove(os.path.join(folder, 'tot_references.csv'))
		except Exception as e:
			continue

for gamma_folder in gamma_folder_list:
	os.chdir(dir_inicial)
	DEF = os.path.join(dir_inicial, "DEFINITIONS.txt")
	tosubmit = os.path.join(dir_inicial, "tosubmit.sh")
	os.chdir(os.path.join(dir_inicial,f"gamma_{gamma_folder}"))

	dir_fuente = os.getcwd()
	folders = []
	folders = ["binary", "binary_ref", "binary_ref/part2"]
	for cell in cell_part:
		folders.append(f"part2/{cell}")
		folders.append(f"part2/{cell}_ref")

	for folder in folders:
		os.makedirs(folder, exist_ok=True)
		if os.path.exists(folder):
			try:
				#shutil.rmtree(folder)
				os.remove(os.path.join(folder, 'tot_references.csv'))
			except Exception as e:
				continue

	lines = read_DEF(DEF)
	gamma = float(gamma_folder.replace('_','.'))
	lines = update_particle_sizes(lines, gamma, R1_np, n_k_bin["part1"], n_k_bin["part2"])
	output_DEF = os.path.join("binary", "DEFINITIONS.txt")
	write_DEF(output_DEF, lines)

	os.chdir("binary")
	DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
	R1_np, R2_np = extract_R_bin(DEF)
	update_cdiva("DEFINITIONS.txt", name_bin, gamma_calc(DEF))
	aL = float(run_command(f'python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {R2_np} {gamma_calc(DEF)} {chain_lenght["part1"]} {chain_lenght["part2"]} {cov["part1"]} {cov["part2"]}'))
	delta_dim_bin = [entry for entry in gamm_delta_dim if entry["gamma"] == gamma]
	process_principal_binario(DEF, name_bin, delta_dim_bin, aL, n_k_bin, tosubmit, dir_fuente, k_aL, gamma, dir_script)
	os.chdir(dir_fuente)

	DEF_part = {}
	for cell in cell_part:
		if flag_reflexion_part == True:
			DEF_part[cell] = os.path.join(dir_script,"references", f"DEFINITIONS_{cell}_reflex.txt")
		else:
			DEF_part[cell] = os.path.join(dir_script,"references", f"DEFINITIONS_{cell}.txt")
	R_part = {"part1": R1_np, "part2": R2_np}
	if os.path.exists(os.path.join(dir_fuente,"part2")):
		dir_fuente = {"part1": os.path.join(dir_inicial,"sim_part1"),"part2": os.path.join(os.getcwd(),"part2")}
		for label in ["part1", "part2"]:
			os.chdir(dir_fuente[label])
			for label_struc in cell_part:
				os.chdir(os.path.join(dir_fuente[label],label_struc))
				dir_fuente_part = dir_fuente[label]
				shutil.copy(DEF_part[label_struc], os.path.join(os.getcwd(),"DEFINITIONS.txt"))
				DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
				lines = read_DEF(DEF)
				size_index = None
				for i, line in enumerate(lines):
					if line.strip() == "!seed":
						size_index = i + 1
						lines[size_index] = seed
						lines[size_index + 1] = seed_lig
					elif line.strip() == "!particle semiaxis x y z in nm":
						size_index = i + 1
						for n in np.arange(0,k_part[label_struc]):
							lines[size_index + n] = f"{R_part[label]} {R_part[label]} {R_part[label]}\n"
					elif line.strip() == "! coverage":
						size_index = i + 1
						for n in np.arange(0,k_part[label_struc]):
							lines[size_index+n] = f"{cov[label]}\n"
					elif line.strip() == "!chains lenght":
						size_index = i + 1
						for n in np.arange(0,k_part[label_struc]):
							lines[size_index+n] = f"{chain_lenght[label]}\n"

				write_DEF(DEF, lines) 

				DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
				R_np = extract_R_part(DEF)
				aL = float(run_command(f'python3 {dir_script}/references/aL_min_{label_struc}.py {R_np} {chain_lenght[label]} {cov[label]}'))
				k_aL_part = 1
				if flag_reflexion_part == True:
					k_aL_part = 2
				delta_dim_part = [entry for entry in part_delta_dim if (entry["part"] == label and entry["cell"] == label_struc)]
				process_principal_part(DEF, delta_dim_part, aL, tosubmit, dir_fuente[label], k_aL_part, dir_script, chain_lenght[label])
					
folder_ref = [os.path.join(dir_inicial,"sim_part1/binary_ref")]
for cell in cell_part:
	folder_ref.append(os.path.join(dir_inicial,f"sim_part1/{cell}_ref"))

for folder in folder_ref:
	os.chdir(folder)
	df = pd.read_csv('tot_references.csv')
	df_fixed = df.drop_duplicates()
	df_fixed.to_csv('tot_references.csv', index=False)

for label in ["part1", "part2"]:
	for gamma_folder in gamma_folder_list:
		dir_fuente = {"part1": os.path.join(dir_inicial,"sim_part1"),"part2": os.path.join(dir_inicial,f'gamma_{gamma_folder}')}
		os.chdir(dir_fuente[label])
		os.chdir("binary_ref")
		references_delta = extract_references("tot_references.csv")
		process_terciario_binario(os.getcwd(), name_bin, references_delta[1:], tosubmit, dir_fuente[label], n_k_bin)
		print(f'created binary {label} gamma = {gamma_folder.replace("_",".")}')

		dir_fuente["part2"] = os.path.join(dir_inicial,f"gamma_{gamma_folder}","part2")
		if os.path.exists(dir_fuente["part2"]):
			for label_struc in cell_part:
				os.chdir(dir_fuente[label])
				DEF = os.path.join(dir_fuente[label], label_struc, "DEFINITIONS.txt")
				os.chdir(os.path.join(dir_fuente[label], f"{label_struc}_ref"))
				references_delta = extract_references("tot_references.csv")
				process_terciario_part(os.getcwd(), references_delta[1:], DEF, tosubmit)

			print(f'created {label} gamma = {gamma_folder.replace("_",".")}')