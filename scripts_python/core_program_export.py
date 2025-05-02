import os
import sys
import subprocess
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import CubicSpline
from collections import defaultdict
from function import run_command, extract_R_bin, extract_cdiva, gamma_calc
from function import path_carpeta, extract_params_init, read_DEF, join_F_csv, join_F_csv_ref
from function_part import extract_R_part
from export_output import process_principal, process_reference_bin, process_principal_part, process_reference_part
from calculate_energy import estimate_part_F, estimate_part_contrib, estimate_bin_F, estimate_bin_contrib, delta_energy_F, delta_energy_US, delta_energy_contrib
from calculate_energy import estimate_bin_volume_overlap, estimate_part_volume_overlap, packing_frac_part, packing_frac_bin

################### INICIO ##################

dir_origin= os.getcwd()
dir_script = os.path.expanduser("~/develop/branch/crystalCF/scripts_python")

params_init = extract_params_init('init_params.txt', False)
name_bin = params_init['name']

flag_reflexion = params_init["flag reflexion binary"]
flag_reflexion_part = params_init["flag reflexion part"]

if flag_reflexion == True:
	if name_bin == "NaCl" or name_bin == "CsCl":
		lines = read_DEF('init_params.txt')
		i = 0
		data = {"n1": None, "n2": None, "R1": None, "num cell bin": None}

		while i < len(lines):
			line = lines[i].strip()
			if line.startswith(("n1", "n2")):
				key, value = line.split()
				data[key] = int(value)
			elif line == "!num cell bin":
				data["num cell bin"] = int(lines[i+1].split()[1])
				i += 1
			i += 1
		n1 = data["n1"]; n2 = data["n2"]; k_bin = data['num cell bin']
	else:
		n1 = params_init['n1']; n2 = params_init['n2']
		k_bin = params_init['num cell bin']
else: 
	n1 = params_init['n1']; n2 = params_init['n2']
	k_bin = params_init['num cell bin']

n = {"part1": n1, "part2": n2}
gamma_list = params_init['gamma list']

gamm_delta_dim = params_init['list gamma delta sum dim']
delta_part = params_init["list delta part"]
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]

gen_curves_flag = params_init["flag generate energy vs aL curves"]
if os.path.isdir(f"results_{name_bin}"):
	shutil.rmtree(f"results_{name_bin}")
os.makedirs(f"results_{name_bin}", exist_ok=True)
final_output = os.path.join(dir_origin,f"results_{name_bin}")

factor_aL_part = {"fcc": 2**(-1.0/6.0), "bcc": 1}

k_reflex_bin = {"kx": 1,"ky": 1,"kz": 1}
if flag_reflexion == True:
	DEF = os.path.join(dir_origin, "DEFINITIONS.txt")
	lines = read_DEF(DEF)
	for i, line in enumerate(lines):
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
			for i, k in enumerate(k_reflex_bin):
				if PBC[2*i] == 3:
					k_reflex_bin[k] = 2
				else:
					k_reflex_bin[k] = 1

################## EXPORTACTION #############
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

F_U = ["F_trans","F_trans_sv","F_vdW"]
F_ST = ["F_conf","F_conf_sv","F_mixs", "F_HS"]
F_name = F_U+F_ST+['F_tot_gcanon']+['volume_overlap']

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
				R1_np, R2_np = extract_R_bin(DEF)
				R = {"part1": R1_np, "part2": R2_np}
				gamma = float(gamma_folder.replace('_','.'))
				aL = float(run_command(f"python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {R2_np}"))
				delta_dim_bin = [entry for entry in gamm_delta_dim if entry["gamma"] == gamma]
				process_principal(output_file, name_bin, R, delta_dim_bin, aL, k_reflex_bin, f_name)
				os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))

				output_file = os.path.join(output_folder, f"{name_bin}_references_{f_name}.csv")
				process_reference_bin(output_file, dir_origin, f_name, R, gamma_folder)

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
						k_reflex_part = 1
						if flag_reflexion_part == True:
							k_reflex_part = 2
						process_principal_part(output_file, label_struc, R_np, delta_list, aL, k_reflex_part, f_name)

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

params_init = extract_params_init('init_params.txt', True)
n1 = params_init['n1']; n2 = params_init['n2']
k_bin = params_init['num cell bin']
n = {"part1": n1, "part2": n2}

for gamma_folder in gamma_folder_list:
	os.chdir(dir_origin)
	os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))
	dir_fuente = os.getcwd()
	DEF = os.path.join(dir_fuente,'binary','DEFINITIONS.txt')
	lines = read_DEF(DEF)
	R1_np, R2_np = extract_R_bin(DEF)
	cdiva_bin = extract_cdiva(DEF)
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

	gamma =  gamma_calc(DEF)
	R = {"part1": R1_np, "part2": R2_np}
	os.chdir("data_analysis")
	df_delta = pd.DataFrame()
	df_contrib = pd.DataFrame()

	keys_list = ["#part", "cell"]+F_U+F_ST
	dict_contrib = {key: [] for key in keys_list}
	dict_delta = {"#part": [],"cell": [], "aL_min": [], "eta": [], "ΔU_min": [], "-ΔST_min": [], "ΔF_min": []}
	V_part_sum = {"part1":{"fcc": [], "bcc": []}, "part2": {"fcc": [], "bcc": []}}
	V_part = {"part1":{"fcc": [], "bcc": []}, "part2": {"fcc": [], "bcc": []}}

	k_reflex_part = 1
	if flag_reflexion_part == True:
		k_reflex_part = 2
	for part in ["part1", "part2"]:
		result_cell = []
		fig, (axp1,axp2) = plt.subplots(ncols = 2, nrows = 1)
		for i, cell in enumerate(cell_part):
			result_cell = estimate_part_F(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, gen_curves_flag, k_reflex_part)
			aL_min = result_cell[0]
			F_part = result_cell[1]
			aL_array = result_cell[2]
			U = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, F_U, aL_array, aL_min, gen_curves_flag, k_reflex_part)
			S = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], n1, n2, F_ST, aL_array, aL_min, gen_curves_flag, k_reflex_part)
			
			V_cell_min, V_sum_cell_min = estimate_part_volume_overlap(part, cell, factor_aL_part[cell], k_part[cell], aL_array, aL_min, k_reflex_part, gen_curves_flag, axp1, axp2)
			V_part_sum[part][cell] =+ V_sum_cell_min
			V_part[part][cell] = V_cell_min
			eta_part = packing_frac_part(part, cell, factor_aL_part[cell], k_part[cell], aL_array, aL_min, k_reflex_part, gen_curves_flag)
			for i, key in enumerate(dict_delta):
				list_params = [part,cell,aL_min, eta_part,U,S,F_part]
				dict_delta[key].append(list_params[i])
		fig.savefig(f"Vol_overlap_{part}.png", format="png")
		plt.close(fig)

	factor_aL_bin = cdiva_bin**(-1.0/3.0)
	result_bin = estimate_bin_F(name_bin, factor_aL_bin, k_bin, n1, n2, ax1, np.round(gamma,2), gen_curves_flag, k_reflex_bin)
	aL_min = result_bin[0]
	F_bin = result_bin[1]
	aL_array = result_bin[2]

	U_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_U, aL_array, aL_min, ax2, np.round(gamma,2), gen_curves_flag, k_reflex_bin)
	S_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_ST, aL_array, aL_min, ax3, np.round(gamma,2), gen_curves_flag, k_reflex_bin)
	V_bin_min, V_sum_bin_min = estimate_bin_volume_overlap(name_bin, factor_aL_bin, k_bin, k_reflex_bin, n1, n2, aL_array, aL_min, gamma, gen_curves_flag)
	eta_bin = packing_frac_bin(name_bin, factor_aL_bin, aL_array, aL_min, k_reflex_bin, gen_curves_flag)
	for i,key in enumerate(dict_delta):
		list_params = ["R2 [nm]", R2_np, "", "", "", "", ""]
		dict_delta[key].append(list_params[i])

	for i,key in enumerate(dict_delta):
		list_params = ["binary",name_bin, aL_min, eta_bin ,U_bin,S_bin,F_bin]
		dict_delta[key].append(list_params[i])

	for key in dict_delta:
		dict_delta[key].append("")

	result = delta_energy_F(dict_delta, cell_part, n1, n2)
	DF = result[0]
	DU = delta_energy_US(dict_delta, "ΔU", cell_min=result[1], n1=n1, n2=n2)
	DS = delta_energy_US(dict_delta, "-ΔST", cell_min=result[1], n1=n1, n2=n2)
	
	for i,key in enumerate(dict_delta):
		list_params = ["gamma",gamma, "","Global X",DU,DS,DF]
		dict_delta[key].append(list_params[i])

	df_delta = pd.DataFrame.from_dict(dict_delta)

	dict_vol = {"#part": [],"cell": [], "#sub_part": [], "Volume overlap [nm³]": []}
	for i, key in enumerate(dict_vol):
		for part in ["part1", "part2"]:
			for k, cell in enumerate(cell_part):
				for l in range(len(V_part[part][cell])):
					list_params = [part,cell,l+1,V_part[part][cell][l][0]]
					dict_vol[key].append(list_params[i])
	for key in dict_vol:
		dict_vol[key].append("")
	for i, key in enumerate(dict_vol):
		for part in ["part1", "part2"]:
			for k, cell in enumerate(cell_part):
				list_params = [part,cell,"Tot",V_part_sum[part][cell][0]]
				dict_vol[key].append(list_params[i])
	
	for key in dict_vol:
		dict_vol[key].append("")
	for i, key in enumerate(dict_vol):
		for part in ["part1", "part2"]:
			for l in range(len(V_bin_min[part])):
				list_params = [name_bin, part,l+1,V_bin_min[part][l][0]]
				dict_vol[key].append(list_params[i])

	for key in dict_vol:
		dict_vol[key].append("")
	for i, key in enumerate(dict_vol):
		list_params = ["binary",name_bin,"Tot",V_sum_bin_min[0]]
		dict_vol[key].append(list_params[i])

	df_vol = pd.DataFrame.from_dict(dict_vol)

	with pd.ExcelWriter(f"results_gamma_{gamma_folder}.xlsx") as writer:
		df_delta.to_excel(writer, sheet_name="Deltas", index=False)
		df_vol.to_excel(writer, sheet_name="Volume Overlap", index=False)

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
eta_values = [[],[],[]]

dfs = []; columnas = 2
for gamma_folder in gamma_folder_list:
	df = pd.read_excel(os.path.join(dir_origin,f"gamma_{gamma_folder}","data_analysis", f"results_gamma_{gamma_folder}.xlsx"),sheet_name='Deltas')
	
	for j, name in enumerate(["part1", "part2", "binary", "Global"]):
		if j<3:
			values = df.loc[df["#part"] == name, ["eta","ΔU_min", "-ΔST_min", "ΔF_min"]]
			gamma = df.loc[df["eta"] == "Global X", ["cell"]].iloc[0]
			if name != "binary":
				eta, DU, DS, DF = values.sort_values(by="ΔF_min").iloc[0]
			else:
				eta, DU, DS, DF = values.iloc[0]

			eta_values[j].append(eta)
		elif j == 3:
			values = df.loc[df["eta"] == "Global X", ["cell", "ΔU_min", "-ΔST_min", "ΔF_min"]]
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

Vol_overlap_dict = defaultdict(list)

parts = ["part1", "part2"]
cell_part = ["fcc", "bcc"]

for gamma_folder in gamma_folder_list:
	path = os.path.join(dir_origin, f"gamma_{gamma_folder}", "data_analysis", f"results_gamma_{gamma_folder}.xlsx")
	df2 = pd.read_excel(path, sheet_name='Volume Overlap')
	df2["#sub_part"] = df2["#sub_part"].astype(str)

	for part in parts:
		for cell in cell_part:
			# Detectar subparts únicos para este (part, cell)
			subset_mask = (df2["#part"] == part) & (df2["cell"] == cell)
			subparts = df2.loc[subset_mask, "#sub_part"].unique()
			for sub in subparts:
				mask = subset_mask & (df2["#sub_part"] == sub)
				values = df2.loc[mask, "Volume overlap [nm³]"].dropna().tolist()
				if values:
					Vol_overlap_dict[(part, cell, sub)].append(values[0])
				else:
					Vol_overlap_dict[(part, cell, sub)].append(float("nan"))
	
	for part in parts:
		subset_mask = (df2["#part"] == name_bin) & (df2["cell"] == part)
		subparts = df2.loc[subset_mask, "#sub_part"].unique()
		for sub in subparts:
			mask = subset_mask & (df2["#sub_part"] == sub)
			values = df2.loc[mask, "Volume overlap [nm³]"].tolist()
			if values:
				Vol_overlap_dict[(name_bin, part, sub)].append(values[0])
			else:
				Vol_overlap_dict[(name_bin, part, sub)].append(float("nan"))

	mask = (df2["#part"] == "binary") & (df2["#sub_part"] == 'Tot')
	values = df2.loc[mask, "Volume overlap [nm³]"].tolist()
	if values:
		Vol_overlap_dict[("binary", name_bin, 'Tot')].append(values[0])
	else:
		Vol_overlap_dict[("binary", name_bin, 'Tot')].append(float("nan"))

	os.chdir(dir_origin)

#datos MD:
ref_MD = pd.read_excel(os.path.join(dir_script,"references","ref_MD_backup.xlsx"), engine="openpyxl")
ref_MD = ref_MD.loc[ref_MD.iloc[:, 0] == name_bin]
gamma_MD = ref_MD[ref_MD.columns[1]]
F_MD = ref_MD[ref_MD.columns[2:]]
y_label = [r'$\Delta$U (k$_{\text{B}}$T)',r'$-T\Delta$S (k$_{\text{B}}$T)',r'$\Delta$F (k$_{\text{B}}$T)']

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
	plt.close()

#####################################################
df_export = pd.DataFrame.from_dict(Vol_overlap_dict, orient="index")
df_export.columns = [f"gamma_{g}" for g in gamma_folder_list]
df_export.reset_index(inplace=True)

df_export.drop(columns=['#part', 'cell', '#sub_part'], errors='ignore', inplace=True)
df_export[['#part', 'cell', '#sub_part']] = pd.DataFrame(df_export['index'].to_list(), index=df_export.index)
df_export.drop(columns=['index'], inplace=True)

df_export = df_export[['#part', 'cell', '#sub_part'] + [c for c in df_export.columns if c not in ['#part', 'cell', '#sub_part']]]

sub_part_sort_map = {'1': 0, '2': 1, 'Tot': 2}
df_export['sub_part_sort'] = df_export['#sub_part'].astype(str).map(sub_part_sort_map)
df_export.sort_values(by=['#part', 'cell', 'sub_part_sort'], inplace=True)
df_export.sort_values(by=['#part','cell'])
df_export.drop(columns=['sub_part_sort'], inplace=True)

output_dict_path = os.path.join(dir_origin, final_output, f"{name_bin}_Vol_overlap_dict.xlsx")
df_export.to_excel(output_dict_path, index=False)

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
		plt.ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=18)
		plt.xlabel(r'$\gamma$',fontsize=18)
		plt.title(f"{name_bin}",fontsize=22,weight='bold')
		plt.legend(fontsize=13)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔU.png',format='png',dpi=300,bbox_inches='tight')
	plt.close()

	DS = pd.DataFrame(DS_values).to_numpy()
	plt.figure(figsize=(8, 6))
	for i in range(len(label)):
		plt.plot(gamma_value,DS[i],marker=marker[i],ms=7, label=label[i])

		plt.axhline(0,ls='--',c='darkgray',zorder=-3)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
		plt.ylabel(r'$-T\Delta$S (k$_{\text{B}}$T)',fontsize=18)
		plt.xlabel(r'$\gamma$',fontsize=18)
		plt.title(f"{name_bin}",fontsize=22,weight='bold')
		plt.legend(fontsize=13)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_-TΔS.png',format='png',dpi=300,bbox_inches='tight')
	plt.close()

DF = pd.DataFrame(DF_values).to_numpy()
plt.figure(figsize=(8, 6))
for i in range(len(label)):
	plt.plot(gamma_value,DF[i],marker=marker[i],ms=7, label=label[i])

plt.axhline(0,ls='--',c='darkgray',zorder=-3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
plt.ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=18)
plt.xlabel(r'$\gamma$',fontsize=18)
plt.title(f"{name_bin}",fontsize=22,weight='bold')
plt.legend(fontsize=13)
plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔF.png',format='png',dpi=300,bbox_inches='tight')
plt.close()

#
dVol = df_export.melt(
    id_vars=['#part', 'cell', '#sub_part'],
    var_name='gamma_label',
    value_name='volume'
)
marker = {"part1": "o", "part2": "v", f"{name_bin}": "d"}
color = {"part1": "tab:blue", "part2": "tab:orange", f"{name_bin}": "tab:green"}
dVol['gamma'] = dVol['gamma_label'].str.replace('gamma_', '', regex=False).str.replace('_', '.').astype(float)
dVol.sort_values(by='gamma', inplace=True)
grouped = dVol.groupby(['#part', 'cell', '#sub_part'],observed=True)

used_labels_ax1 = set(); used_labels_ax2 = set()
fig, (ax1,ax2) = plt.subplots(ncols = 2,nrows = 1, figsize=(12,6),constrained_layout=True)
for (part, cell, sub), group in grouped:
	if sub != 'Tot':
		label = part if part not in used_labels_ax1 else None
		ax1.plot(gamma_value,group['volume'],c=color[part],marker=marker[part], label=label)
		used_labels_ax1.add(part)
	else:
		if part == "binary":
			part = name_bin
		label = part if part not in used_labels_ax2 else None
		ax2.plot(gamma_value,group['volume'],c=color[part],marker=marker[part], label=label)
		used_labels_ax2.add(part)
	
for ax in [ax1,ax2]:
	ax.tick_params(axis='both', labelsize=14)
	ax.set_xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
	ax.set_ylabel(r'Vol. overlap [nm$^3$]',fontsize=16)
	ax.set_xlabel(r'$\gamma$',fontsize=18)
	ax.legend(fontsize=13)

plt.savefig(f'{final_output}/{name_bin}_contrib_volume_overlap.png',format='png',dpi=300,bbox_inches='tight')

# eta
label = ["part1", "part2", "binary", "global"]
marker = ["o", "v", "d", "s"]
color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

DF = pd.DataFrame(DF_values).to_numpy()
eta = pd.DataFrame(eta_values).to_numpy()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

# Primer subplot: ΔF
for i in range(len(label)):
    ax1.plot(gamma_value, DF[i], color=color[i], marker=marker[i], ms=7, label=label[i])
ax1.axhline(0, ls='--', c='darkgray', zorder=-3)
ax1.set_xlim(min(gamma_value) - 0.05, max(gamma_value) + 0.05)
ax1.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)', fontsize=18)
ax1.tick_params(axis='both', labelsize=14)
ax1.set_title(f"{name_bin}", fontsize=22, weight='bold')
ax1.legend(fontsize=13, loc='best')

# Segundo subplot: η
for i in range(len(label)-1):  # Asumes que eta tiene len(label)-1 filas
    ax2.plot(gamma_value, eta[i], ls='--', marker=marker[i], color=color[i], label=label[i], linewidth=2)
ax2.set_xlim(min(gamma_value) - 0.05, max(gamma_value) + 0.05)
ax2.set_xlabel(r'$\gamma$', fontsize=18)
ax2.set_ylabel(r'$\eta$', fontsize=18)
ax2.tick_params(axis='both', labelsize=14)
ax2.legend(fontsize=13, loc='best')

plt.tight_layout()
plt.savefig(f'{final_output}/{name_bin}_DF_and_eta_results.png', format='png', dpi=300, bbox_inches='tight')
plt.close()

F_plot = DF_values[3]
eta_values = np.array(eta_values)
eta_plot = eta_values[2]

plt.figure(figsize=(8, 6))
plt.plot(eta_plot,F_plot,ls='none',marker='s',color='red',ms=7,label='MoltCF')
plt.axhline(0,ls='--',c='darkgray',zorder=-3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=18)
plt.xlabel(rf'$\eta\ {name_bin} $',fontsize=18)
plt.title(f"{name_bin}",fontsize=22,weight='bold')
plt.legend(fontsize=13)
plt.savefig(f'{final_output}/{name_bin}_results_F_vs_eta.png',format='png',dpi=300,bbox_inches='tight')
plt.close()
