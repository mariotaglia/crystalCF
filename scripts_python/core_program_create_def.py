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
from function import extract_R_bin, extract_references, update_particle_sizes, update_R1, extract_definitions, update_nseg_cov
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

Nseg_list = params_init['nseg list']
Cov_list = params_init['coverage list']
Nseg_cov_delta_dim = params_init['list nseg coverage delta sum dim']
k_bin = params_init['num cell bin']
n_k_bin = {"part1": n1*k_bin, "part2": n2*k_bin}
flag_reflexion = params_init["flag reflexion binary"]
flag_reflexion_part = params_init["flag reflexion part"]

delta_part = params_init["list delta part"]
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]
nseg_folder_list = ["{:.0f}".format(n).replace('.','_') for n in Nseg_list]
cov_folder_list = ["{:.2f}".format(c).replace('.','_') for c in Cov_list]

update_R1("DEFINITIONS.txt", n_k_bin["part1"], R1_np)
DEF = os.path.join(dir_inicial, "DEFINITIONS.txt")
lines = read_DEF(DEF)
for i, line in enumerate(lines):
	if line.strip() == "!seed":
		size_index = i + 1
		seed = lines[size_index]
		seed_lig = lines[size_index + 1]
		break


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

for nseg_folder in nseg_folder_list:
	os.makedirs(f'nseg_{nseg_folder}', exist_ok=True)
	os.chdir(f'nseg_{nseg_folder}')

	for cov_folder in cov_folder_list:
		os.makedirs(f'cov_{cov_folder}', exist_ok=True)
		os.chdir(f'cov_{cov_folder}')
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
		os.chdir(os.path.join(dir_inicial,f'nseg_{nseg_folder}'))

	os.chdir(dir_inicial)

for nseg_folder in nseg_folder_list:
	os.chdir(dir_inicial)
	for cov_folder in cov_folder_list:
		DEF = os.path.join(dir_inicial, "DEFINITIONS.txt")
		tosubmit = os.path.join(dir_inicial, "tosubmit.sh")
		os.chdir(os.path.join(dir_inicial,f"nseg_{nseg_folder}",f"cov_{cov_folder}"))

		dir_fuente = os.getcwd()
		folders = ["binary"]

		for folder in folders:
			os.makedirs(folder, exist_ok=True)
			if os.path.exists(folder):
				try:
					#shutil.rmtree(folder)
					os.remove(os.path.join(folder, 'tot_references.csv'))
				except Exception as e:
					continue

		lines = read_DEF(DEF)
		nseg = float(nseg_folder.replace('_','.'))
		cov = float(cov_folder.replace('_','.'))
		lines = update_nseg_cov(lines, nseg, cov, n_k_bin["part1"], n_k_bin["part2"])
		output_DEF = os.path.join("binary", "DEFINITIONS.txt")
		write_DEF(output_DEF, lines)

		os.chdir("binary")
		DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")

		update_cdiva("DEFINITIONS.txt", name_bin, 1, flag_reflexion)
		aL = float(run_command(f'python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {cov} {nseg}'))
		delta_dim_bin = [entry for entry in Nseg_cov_delta_dim if entry["nseg"] == nseg and entry["coverage"] == cov]
		process_principal_binario(DEF, name_bin, delta_dim_bin, aL, n_k_bin, tosubmit, dir_fuente, k_aL, nseg, cov)
		os.chdir(dir_fuente)

		DEF_part = {}
		for cell in cell_part:
			if flag_reflexion_part == True:
				DEF_part[cell] = os.path.join(dir_script,"references", f"DEFINITIONS_{cell}_reflex.txt")
			else:
				DEF_part[cell] = os.path.join(dir_script,"references", f"DEFINITIONS_{cell}.txt")
		
		R_part = {"part1": R1_np}
		
		dir_fuente_part = {"part1": os.path.join(dir_fuente,"sim_part1")}
		os.chdir(dir_fuente_part["part1"])
		for label_struc in cell_part:
			os.chdir(os.path.join(dir_fuente_part["part1"],label_struc))
			shutil.copy(DEF_part[label_struc], os.path.join(os.getcwd(),"DEFINITIONS.txt"))
			DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
			lines = read_DEF(DEF)
			size_index = None
			for i, line in enumerate(lines):
				if line.strip() == "!seed":
					size_index = i + 1
					lines[size_index] = seed
					lines[size_index + 1] = seed_lig
				elif line.strip() == "!properties of ligand chains":
					size_index = i + 1
					lines[size_index] = f"long {str(int(nseg))}\n"
				elif line.strip() == "!particle semiaxis x y z in nm":
					size_index = i + 1
					for n in np.arange(0,k_part[label_struc]):
						lines[size_index + n] = f"{R_part["part1"]} {R_part["part1"]} {R_part["part1"]}\n"
				elif line.strip() == "! coverage":
					size_index = i + 1
					for n in np.arange(0,k_part[label_struc]):
						lines[size_index + n] = f"{cov}\n"
					break

			write_DEF(DEF, lines) 

			DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
			R_np = extract_R_part(DEF)
			aL = float(run_command(f'python3 {dir_script}/references/aL_min_{label_struc}.py {R_np} {nseg} {cov}'))
			k_aL_part = 1
			if flag_reflexion_part == True:
				k_aL_part = 2
			process_principal_part(DEF, delta_part[label_struc], aL, tosubmit, dir_fuente_part["part1"], k_aL_part)

		folder_ref = [os.path.join(dir_fuente,"sim_part1/binary_ref")]
		for cell in cell_part:
			folder_ref.append(os.path.join(dir_fuente,f"sim_part1/{cell}_ref"))

		for folder in folder_ref:
			os.chdir(folder)
			df = pd.read_csv('tot_references.csv')
			df_fixed = df.drop_duplicates()
			df_fixed.to_csv('tot_references.csv', index=False)

		for label in ["part1"]:
			dir_fuente_part = {"part1": os.path.join(dir_fuente,"sim_part1")}
			os.chdir(dir_fuente_part[label])
			os.chdir("binary_ref")
			references_delta = extract_references("tot_references.csv")
			process_terciario_binario(os.getcwd(), name_bin, references_delta[1:], tosubmit, dir_fuente_part[label], n_k_bin)
			print(f'created {name_bin} nseg = {nseg_folder.replace("_",".")} cov = {cov_folder.replace("_",".")}')

			for label_struc in cell_part:
				os.chdir(dir_fuente_part[label])
				DEF = os.path.join(dir_fuente_part[label], label_struc, "DEFINITIONS.txt")
				os.chdir(os.path.join(dir_fuente_part[label], f"{label_struc}_ref"))
				references_delta = extract_references("tot_references.csv")
				process_terciario_part(os.getcwd(), references_delta[1:], DEF, tosubmit)

			print(f'created {label} nseg = {nseg_folder.replace("_",".")} cov = {cov_folder.replace("_",".")}')
