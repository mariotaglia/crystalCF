import os
import subprocess
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import CubicSpline
from collections import defaultdict
from function import run_command, extract_R_bin, gamma_calc
from function import path_carpeta, extract_params_init, read_DEF, join_F_csv, join_F_csv_ref
from function_part import extract_R_part
from export_output import process_principal, process_reference_bin, process_principal_part, process_reference_part
from calculate_energy import estimate_part_F, estimate_part_contrib, estimate_bin_F, estimate_bin_contrib, delta_energy_F, delta_energy_US, delta_energy_contrib

################### INICIO ##################

dir_origin= os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python")

params_init = extract_params_init('init_params.txt')
name_bin = params_init['name']
n1 = params_init['n1']; n2 = params_init['n2']
n = {"part1": n1, "part2": n2}

gamma_list = params_init['gamma list']
k_bin = params_init['num cell bin']
gamm_delta_dim = params_init['list gamma delta sum dim']
factor_aL_bin = params_init['aL cell bin factor']

delta_part = params_init["list delta part"]
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]
factor_aL_part = params_init['aL cell part factor']

gen_curves_flag = params_init["flag generate energy vs aL curves"]
if os.path.isdir(f"results_{name_bin}"):
	shutil.rmtree(f"results_{name_bin}")
os.makedirs(f"results_{name_bin}", exist_ok=True)
final_output = os.path.join(dir_origin,f"results_{name_bin}")

################## EXPORTACTION #############
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

F_U = ["F_trans","F_trans_sv","F_vdW"]
F_ST = ["F_conf","F_conf_sv","F_mixs", "F_HS"]
F_name = F_U+F_ST+['F_tot_gcanon']

while True:
	check_extract = input("¿Quiere extraer los F.data? (s/n): ").strip().lower()
	if check_extract in ["s", "si"]:
		print("Extrayendo F.data...")
		for gamma_folder in gamma_folder_list:
			os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))

			dir_fuente = os.getcwd()
			if os.path.isdir("data_analysis"):
				shutil.rmtree("data_analysis")

			os.makedirs("data_analysis", exist_ok=True)
			output_folder = os.path.join(dir_fuente,"data_analysis")

			for f_name in F_name:
				os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))
				os.chdir("binary")
				output_file = os.path.join(output_folder, f"{name_bin}_results_{f_name}.csv")
				DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
				R1_np, R2_np = extract_R_bin(DEF, n1*k_bin)
				
				gamma = float(gamma_folder.replace('_','.'))
				aL = float(run_command(f"python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {R2_np}"))
				delta_dim_bin = [entry for entry in gamm_delta_dim if entry["gamma"] == gamma]
				process_principal(output_file, delta_dim_bin, aL, f_name)
				os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))

				output_file = os.path.join(output_folder, f"{name_bin}_references_{f_name}.csv")
				process_reference_bin(output_file, dir_origin, f_name)

				dir_fuente = {"part1": os.path.join(dir_origin,"sim_part1"),"part2": os.path.join(os.getcwd(),"part2")}
				for label in ["part1", "part2"]:
					for label_struc in cell_part:
						os.chdir(dir_fuente[label])
						output_file = os.path.join(output_folder, f"{label}_results_{f_name}.csv")
						DEF = os.path.join(dir_fuente[label], label_struc, "DEFINITIONS.txt")
						os.chdir(f"{label_struc}")
						R_np = extract_R_part(DEF)
						aL = float(run_command(f"python3 {dir_script}/references/aL_min_{label_struc}.py {R_np}"))
						delta_list = delta_part[label_struc]
						process_principal_part(output_file, label_struc, delta_list, aL, f_name)

						os.chdir(os.path.join(dir_fuente[label], f"{label_struc}_ref"))
						output_file = os.path.join(output_folder, f"{label}_references_{f_name}.csv")
						base_folder = os.path.join(dir_fuente[label], f"{label_struc}_ref")
						process_reference_part(output_file, base_folder, cell_part, label_struc, f_name)

			join_F_csv(output_folder, name_bin, True)
			join_F_csv_ref(output_folder, name_bin)
			for label in ["part1","part2"]:
				join_F_csv(output_folder, label, False)
				join_F_csv_ref(output_folder, label)

			print(f"Exportado gamma {gamma_folder.replace('_','.')}")
		break

	elif check_extract in ["n", "no"]:
		print("Salteado la exportacion.")
		break

	else:
		print("Respuesta no válida. Por favor, ingrese 's' o 'n'.")

##################### ESTIMACIONES ###############################
os.chdir(dir_origin)

import matplotlib.pyplot as plt
fig1, ax1 = plt.subplots(); fig2, ax2 = plt.subplots(); fig3, ax3 = plt.subplots()

for gamma_folder in gamma_folder_list:
	os.chdir(dir_origin)
	os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))
	dir_fuente = os.getcwd()
	DEF = os.path.join(dir_fuente,'binary','DEFINITIONS.txt')
	lines = read_DEF(DEF)
	size_index = None
	for i, line in enumerate(lines):
		if line.strip() == "!properties of ligand chains":
			size_index = i+1
			nseg = float(lines[size_index].split()[1])
			size_index = None
		if line.startswith("! coverage"):
			size_index = i+1
			sigma = float(lines[size_index].split()[0])
		if line.startswith("!volume"):
			size_index = i+1
			vsol = float(lines[size_index].split()[1])
		if line.strip() == "!particle semiaxis x y z in nm":
			size_index = i + 1
			R1_np = float(lines[size_index].split()[0])  # Tomar el primer valor de la línea
			R2_np = float(lines[size_index + n1].split()[0])  # Tomar el primer valor de la línea 
			size_index = None

	gamma =  gamma_calc(DEF, n1*k_bin)
	R = {"part1": R1_np, "part2": R2_np}
	os.chdir("data_analysis")
	df_delta = pd.DataFrame()
	df_contrib = pd.DataFrame()

	keys_list = ["#part", "cell"]+F_U+F_ST
	dict_contrib = {key: [] for key in keys_list}
	dict_delta = {"#part": [],"cell": [], "aL_min": [], "ΔU_min": [], "-ΔST_min": [], "ΔF_min": []}

	for part in ["part1", "part2"]:
		result_cell = []
		for i, cell in enumerate(cell_part):
			result_cell = estimate_part_F(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, gen_curves_flag)
			aL_min = result_cell[0]
			F_part = result_cell[1]
			aL_array = result_cell[2]
			U = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, F_U, aL_array, aL_min, gen_curves_flag)
			S = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, F_ST, aL_array, aL_min, gen_curves_flag)

			dict_delta["#part"].append(part)
			dict_delta["cell"].append(cell)
			dict_delta["aL_min"].append(aL_min)
			dict_delta["ΔU_min"].append(U)
			dict_delta["-ΔST_min"].append(S)
			dict_delta["ΔF_min"].append(F_part)

	result_bin = estimate_bin_F(name_bin, factor_aL_bin, k_bin, n1, n2, ax1, np.round(gamma,2), gen_curves_flag)
	aL_min = result_bin[0]
	F_bin = result_bin[1]
	aL_array = result_bin[2]

	U_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_U, aL_array, aL_min, ax2, np.round(gamma,2), gen_curves_flag)
	S_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_ST, aL_array, aL_min, ax3, np.round(gamma,2), gen_curves_flag)

	for key in dict_delta:
		dict_delta[key].append("")

	dict_delta["#part"].append("binary")
	dict_delta["cell"].append(name_bin)
	dict_delta["aL_min"].append(aL_min)
	dict_delta["ΔU_min"].append(U_bin)
	dict_delta["-ΔST_min"].append(S_bin)
	dict_delta["ΔF_min"].append(F_bin)

	for key in dict_delta:
		dict_delta[key].append("")

	result = delta_energy_F(dict_delta, cell_part, n1, n2)
	DF = result[0]
	DU = delta_energy_US(dict_delta, "ΔU", cell_min=result[1], n1=n1, n2=n2)
	DS = delta_energy_US(dict_delta, "-ΔST", cell_min=result[1], n1=n1, n2=n2)
	dict_delta["#part"].append("gamma")
	dict_delta["cell"].append(gamma)
	dict_delta["aL_min"].append("Global")
	dict_delta["ΔU_min"].append(DU)
	dict_delta["-ΔST_min"].append(DS)
	dict_delta["ΔF_min"].append(DF)

	df_delta = pd.DataFrame.from_dict(dict_delta)

	with pd.ExcelWriter(f"results_gamma_{gamma_folder}.xlsx") as writer:
		df_delta.to_excel(writer, sheet_name="Deltas", index=False)

if gen_curves_flag == True:
	for ax in [ax1,ax2,ax3]:
		ax.set_xlabel(r'a$_{\text{L}}$',fontsize=16)
		ax.legend()
	for fig in [fig1,fig2,fig3]:
		fig.suptitle(f'{name_bin} binary',fontsize=18)
	ax1.set_ylabel(r'$\Delta$F (k$_{\text{b}}$T)',fontsize=16)
	ax2.set_ylabel(r'$\Delta$U (k$_{\text{b}}$T)',fontsize=16)
	ax3.set_ylabel(r'-T$\Delta$S (k$_{\text{b}}$T)',fontsize=16)

	fig1.savefig(f"{final_output}/F_binary.png", format="png", dpi=300)
	fig2.savefig(f"{final_output}/U_binary.png", format="png", dpi=300)
	fig3.savefig(f"{final_output}/S_binary.png", format="png", dpi=300)
	for fig in [fig1,fig2,fig3]:
		plt.close(fig)

####################### PLOT ###########################
os.chdir(dir_origin)

gamma_value = []
DU_global = []
DS_global = []
DF_global = []

DU_values = [[],[],[],[]]
DS_values = [[],[],[],[]]
DF_values = [[],[],[],[]]

dfs = []; columnas = 2
for gamma_folder in gamma_folder_list:
	df = pd.read_excel(os.path.join(os.path.join(dir_origin,f"gamma_{gamma_folder}","data_analysis", f"results_gamma_{gamma_folder}.xlsx")))
	
	for j, name in enumerate(["part1", "part2", "binary", "Global"]):
		if j<3:
			values = df.loc[df["#part"] == name, ["ΔU_min", "-ΔST_min", "ΔF_min"]]
			gamma = df.loc[df["aL_min"] == "Global", ["cell"]].iloc[0]
			if name != "binary":
				DU, DS, DF = values.sort_values(by="ΔF_min").iloc[0]
			else:
				DU, DS, DF = values.iloc[0]

		elif j == 3:
			values = df.loc[df["aL_min"] == "Global", ["cell", "ΔU_min", "-ΔST_min", "ΔF_min"]]
			gamma, DU, DS, DF = values.iloc[0]
		
		DU_values[j].append(DU)
		DS_values[j].append(DS)
		DF_values[j].append(DF)

	gamma_value.append(gamma)

	os.chdir(dir_origin)

	df.insert(0, "Gamma", gamma_folder)
	# Agregar fila vacía al final
	df = pd.concat([df, pd.DataFrame([[""] * len(df.columns)], columns=df.columns)], ignore_index=True)
	dfs.append(df)

#datos MD:
ref_MD = pd.read_excel(os.path.join(dir_script,"references","ref_MD_backup.xlsx"), engine="openpyxl")
ref_MD = ref_MD.loc[ref_MD.iloc[:, 0] == name_bin]
gamma_MD = ref_MD[ref_MD.columns[1]]
F_MD = ref_MD[ref_MD.columns[2:]]
y_label = [r'$\Delta$U (k$_{\text{b}}$T)',r'$-T\Delta$S (k$_{\text{b}}$T)',r'$\Delta$F (k$_{\text{b}}$T)']

F_plot = [DU_values[3],DS_values[3],DF_values[3]]
for i, (lista, F) in enumerate(zip(F_plot,["ΔU", "-TΔS", "ΔF"])):
	plt.figure(figsize=(8, 6))
	plt.plot(gamma_value,F_plot[i],ls='none',marker='s',color='red',ms=7,label='MoltCF')
	plt.scatter(gamma_MD,F_MD[F],marker='o',color='purple',s=50,label='MD (MD)',zorder=10)

	plt.axhline(0,ls='--',c='darkgray',zorder=-3)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xlim(min(gamma_list)-0.05,max(gamma_list)+0.05)
	plt.ylabel(y_label[i],fontsize=18)
	plt.xlabel(r'$\gamma$',fontsize=18)
	plt.title(f"{name_bin}",fontsize=22,weight='bold')
	plt.legend(fontsize=13)
	plt.savefig(f'{final_output}/{name_bin}_results_{F}.png',format='png',dpi=300,bbox_inches='tight')

#####################################################

if len(dfs) % columnas != 0:
    empty_df = pd.DataFrame([[""] * len(dfs[0].columns)], columns=dfs[0].columns)
    dfs.append(empty_df)
max_filas = max(df.shape[0] for df in dfs)

for i in range(len(dfs)):
    filas_actuales = dfs[i].shape[0]
    if filas_actuales < max_filas:
        filas_faltantes = max_filas - filas_actuales
        filas_vacias = pd.DataFrame([[""] * len(dfs[i].columns)] * filas_faltantes, columns=dfs[i].columns)
        dfs[i] = pd.concat([dfs[i], filas_vacias], ignore_index=True)

tabla_final = []

for i in range(0, len(dfs), columnas):
    fila_data = pd.concat(dfs[i:i + columnas], axis=1)
    fila_data.insert(len(dfs[0].columns), "", "")
    tabla_final.append(fila_data)
    
tabla_final_df = pd.concat(tabla_final, ignore_index=True)

output_path = os.path.join(dir_origin, f"{final_output}/{name_bin}_table_results.xlsx")
tabla_final_df.to_excel(output_path, index=False)

########## graficos U, S ##########
label = ["part1", "part2", "binary", "global"]
marker = ["o", "v", "d", "s"]
if gen_curves_flag == True:
	DU = pd.DataFrame(DU_values).to_numpy()
	plt.figure(figsize=(8, 6))
	for i in range(len(label)):
		plt.plot(gamma_value,DU[i],marker=marker[i],ms=7, label=label[i])

		plt.axhline(0,ls='--',c='darkgray',zorder=-3)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
		plt.ylabel(r'$\Delta$U (k$_{\text{b}}$T)',fontsize=18)
		plt.xlabel(r'$\gamma$',fontsize=18)
		plt.title(f"{name_bin}",fontsize=22,weight='bold')
		plt.legend(fontsize=13)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔU.png',format='png',dpi=300,bbox_inches='tight')

	DS = pd.DataFrame(DS_values).to_numpy()
	plt.figure(figsize=(8, 6))
	for i in range(len(label)):
		plt.plot(gamma_value,DS[i],marker=marker[i],ms=7, label=label[i])

		plt.axhline(0,ls='--',c='darkgray',zorder=-3)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
		plt.ylabel(r'$-T\Delta$S (k$_{\text{b}}$T)',fontsize=18)
		plt.xlabel(r'$\gamma$',fontsize=18)
		plt.title(f"{name_bin}",fontsize=22,weight='bold')
		plt.legend(fontsize=13)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_-TΔS.png',format='png',dpi=300,bbox_inches='tight')

DF = pd.DataFrame(DF_values).to_numpy()
plt.figure(figsize=(8, 6))
for i in range(len(label)):
	plt.plot(gamma_value,DF[i],marker=marker[i],ms=7, label=label[i])

	plt.axhline(0,ls='--',c='darkgray',zorder=-3)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
	plt.ylabel(r'$\Delta$F (k$_{\text{b}}$T)',fontsize=18)
	plt.xlabel(r'$\gamma$',fontsize=18)
	plt.title(f"{name_bin}",fontsize=22,weight='bold')
	plt.legend(fontsize=13)
	plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔF.png',format='png',dpi=300,bbox_inches='tight')