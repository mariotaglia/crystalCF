import os
import re
import sys
import subprocess
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import CubicSpline
from collections import defaultdict
from function import run_command, extract_R_bin, extract_cdiva, gamma_calc, vol_tot_bin
from function import path_carpeta, extract_params_init, read_DEF, join_F_csv, join_F_csv_ref
from function_part import extract_R_part, vol_tot_part
from export_output import process_principal, process_reference_bin, process_principal_part, process_reference_part
from calculate_energy import estimate_part_F,estimate_part_F_pair, estimate_part_contrib, estimate_bin_F, estimate_bin_F_pair, estimate_bin_contrib, delta_energy_F, delta_energy_US, delta_energy_contrib

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='cm10', weight='bold', size=16)
import matplotlib as mpl
mpl.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}"
})

def to_latex_formula(name):
    """
    Convierte una fórmula química tipo texto como 'NaCl', 'AlB2' a LaTeX: NaCl, AlB$_2$
    """
    # Busca los elementos con su posible subíndice (e.g., AlB2 -> [('Al', ''), ('B', '2')])
    match = re.match(r'^([a-z]+)?([A-Z][A-Za-z0-9_]+)$', name)
    if not match:
        return name  # Devuelve tal cual si no coincide el formato esperado

    prefix, formula = match.groups()
    prefix = prefix or ""

    # Reemplazar subíndices tipo _6 o 6 después de letras
    # Ej: AB_6 o AB6 → AB$_6$
    formula = re.sub(r'_(\d+)', r'$_{\1}$', formula)
    formula = re.sub(r'([A-Za-z])(\d+)', r'\1$_{\2}$', formula)

    return f"{prefix} {formula}".strip()

################### INICIO ##################

dir_origin= os.getcwd()
dir_script = os.path.expanduser("~/develop/crystalCF/scripts_python")

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
	
n_k_bin = {"part1": n1*k_bin, "part2": n2*k_bin}
n = {"part1": n1, "part2": n2}
gamma_list = params_init['gamma list']

gamm_delta_dim = params_init['list gamma delta sum dim']
part_delta_dim = params_init['list part delta sum dim']
cell_part = params_init["cell part"]
k_part = params_init["num cell part"]

gen_curves_flag = params_init["flag generate energy vs aL curves"]
if os.path.isdir(f"results_{name_bin}"):
	shutil.rmtree(f"results_{name_bin}")
os.makedirs(f"results_{name_bin}", exist_ok=True)
final_output = os.path.join(dir_origin,f"results_{name_bin}")

cdiva_fcc = np.sqrt(2)
factor_aL_part = {"fcc": cdiva_fcc/np.power(cdiva_fcc,(1./3.)), "bcc": 1}
n1_k_bin = n_k_bin["part1"]; n2_k_bin = n_k_bin["part2"]
k_aL = {"kx": 1,"ky": 1,"kz": 1}
DEF = os.path.join(dir_origin, "DEFINITIONS.txt")
lines = read_DEF(DEF)
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

for i, line in enumerate(lines):
	nseg = []
	for i, line in enumerate(lines):
		if line.strip() == "!properties of ligand chains":
			nseg = float(lines[i+1].split()[1])
			break

for label in ['part1','part2']:
	if chain_lenght[label] == None:
		chain_lenght[label] = nseg

if flag_reflexion == True:
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
			for i, k in enumerate(k_aL):
				if PBC[2*i] == 3:
					k_aL[k] = 2
				else:
					k_aL[k] = 1

################## EXPORTACTION #############
gamma_folder_list = ["{:.3f}".format(g).replace('.','_') for g in gamma_list]

F_U = ["F_trans","F_trans_sv","F_vdW"]
F_ST = ["F_conf","F_conf_sv","F_mixs", "F_HS"]
F_name = F_U+F_ST+['F_tot_gcanon']#+["F_pairwise"]

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
				R = {"part1": np.max([R1_np, R2_np]), "part2": np.min([R1_np, R2_np])}
				R1_np = R["part1"]; R2_np = R["part2"]
				gamma = float(gamma_folder.replace('_','.'))
				aL = float(run_command(f'python3 {dir_script}/references/aL_estimate_bin.py {name_bin} {R1_np} {R2_np} {gamma_calc(DEF)} {chain_lenght["part1"]} {chain_lenght["part2"]} {cov["part1"]} {cov["part2"]}'))
				delta_dim_bin = [entry for entry in gamm_delta_dim if entry["gamma"] == gamma]
				process_principal(output_file, name_bin, R, delta_dim_bin, aL, k_aL, f_name)
				os.chdir(os.path.join(dir_origin,f"gamma_{gamma_folder}"))
				if not f_name == "F_pairwise":
					output_file = os.path.join(output_folder, f"{name_bin}_references_{f_name}.csv")
					if name_bin != "Li3Bi" and name_bin != "NaZn13":
						process_reference_bin(output_file, dir_origin, f_name, R, gamma_folder)
					else:
						if not os.path.isfile(output_file):
							with open(output_file, "w") as out_file:
								out_file.write("#part,radius [nm],delta,dimx,dimy,dimz,F_reference\n")
				dir_fuente = {"part1": os.path.join(dir_origin,"sim_part1"),"part2": os.path.join(os.getcwd(),"part2")}

				for label in ["part1", "part2"]:
					for label_struc in cell_part:
						os.chdir(dir_fuente[label])
						output_file = os.path.join(output_folder, f"{label}_results_{f_name}.csv")
						DEF = os.path.join(dir_fuente[label], label_struc, "DEFINITIONS.txt")
						os.chdir(f"{label_struc}")
						R_np = extract_R_part(DEF)
						aL = float(run_command(f'python3 {dir_script}/references/aL_min_{label_struc}.py {R_np} {chain_lenght[label]} {cov[label]}'))
						
						k_aL_part = 1
						if flag_reflexion_part == True:
							k_aL_part = 2

						delta_dim_part = [entry for entry in part_delta_dim if (entry["part"] == label and entry["cell"] == label_struc)]
						process_principal_part(output_file, label_struc, R_np, delta_dim_part, aL, k_aL_part, f_name)

						os.chdir(os.path.join(dir_fuente[label], f"{label_struc}_ref"))
						output_file = os.path.join(output_folder, f"{label}_references_{f_name}.csv")
						base_folder = os.path.join(dir_fuente[label], f"{label_struc}_ref")
						if not f_name == "F_pairwise":
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
v_pol_part = {"part1": 0.028,"part2": 0.028}
v_pol_bin = {"NaCl": 0.028}
import matplotlib.pyplot as plt
fig1, ax1 = plt.subplots(figsize=(8, 6)); fig2, ax2 = plt.subplots(figsize=(8, 6)); fig3, ax3 = plt.subplots(figsize=(8, 6))
fig4, ax4 = plt.subplots(figsize=(8, 6))
params_init = extract_params_init('init_params.txt', True)
n1 = params_init['n1']; n2 = params_init['n2']
k_bin = params_init['num cell bin']
n = {"part1": n1, "part2": n2}
cell_bin_factor = params_init["cell bin factor"]
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
	dict_delta = {"#part": [],"cell": [], "aL_min": [], "ΔU_min": [], "-ΔST_min": [], "ΔF_min": []}

	k_aL_part = 1
	if flag_reflexion_part == True:
		k_aL_part = 2
	for part in ["part1", "part2"]:
		result_cell = []
		for i, cell in enumerate(cell_part):
			result_cell = estimate_part_F(part, cell, factor_aL_part[cell], n[part], k_part[cell], gen_curves_flag, k_aL_part)
			#result_cell_pairwise = estimate_part_F_pair(part, cell, factor_aL_part[cell], n[part], k_part[cell],vol_tot_part(cell,R[part],chain_lenght[part],cov[part],v_pol_part[part]), gen_curves_flag, k_aL_part, np.round(gamma,2))
			aL_min = result_cell[0]
			F_part = result_cell[1]
			aL_array = result_cell[2]
			#aL_min = result_cell_pairwise[0]
			#F_part = result_cell_pairwise[1]
			#aL_array = result_cell_pairwise[2]
			U = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], F_U, aL_array, aL_min, gen_curves_flag, k_aL_part)
			S = estimate_part_contrib(part, cell, factor_aL_part[cell], n[part], k_part[cell], F_ST, aL_array, aL_min, gen_curves_flag, k_aL_part)
			#U = 0; S = 0
			for i,key in enumerate(dict_delta):
				list = [part,cell,aL_min,U,S,F_part]
				dict_delta[key].append(list[i])

	factor_aL_bin = cell_bin_factor*np.power(cdiva_bin,(-1.0/3.0))
	result_bin = estimate_bin_F(name_bin, factor_aL_bin, k_bin, n1, n2, ax1, np.round(gamma,2), gen_curves_flag, k_aL)
	#result_bin_pair = estimate_bin_F_pair(name_bin, factor_aL_bin, k_bin, n1, n2, vol_tot_bin(name_bin,R1_np,R2_np,chain_lenght["part1"],chain_lenght["part2"],cov["part1"],cov["part2"],0.028), ax4, np.round(gamma,2), gen_curves_flag, k_aL)
	aL_min = result_bin[0]
	F_bin = result_bin[1]
	aL_array = result_bin[2]
	#aL_min = result_bin_pair[0]
	#F_bin = result_bin_pair[1]
	#aL_array = result_bin_pair[2]
	U_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_U, aL_array, aL_min, ax2, np.round(gamma,2), gen_curves_flag, k_aL)
	S_bin = estimate_bin_contrib(name_bin, factor_aL_bin, k_bin, n1, n2, F_ST, aL_array, aL_min, ax3, np.round(gamma,2), gen_curves_flag, k_aL)
	#U_bin = 0; S_bin = 0
	for i,key in enumerate(dict_delta):
		list = ["R2 [nm]", R2_np, "", "", "", ""]
		dict_delta[key].append(list[i])

	for i,key in enumerate(dict_delta):
		list = ["binary",name_bin,aL_min,U_bin,S_bin,F_bin]
		dict_delta[key].append(list[i])

	for key in dict_delta:
		dict_delta[key].append("")

	result = delta_energy_F(dict_delta, cell_part, n1, n2)
	DF = result[0]
	DU = delta_energy_US(dict_delta, "ΔU", cell_min=result[1], n1=n1, n2=n2)
	DS = delta_energy_US(dict_delta, "-ΔST", cell_min=result[1], n1=n1, n2=n2)
	
	for i,key in enumerate(dict_delta):
		list = ["gamma",gamma,"Global X",DU,DS,DF]
		dict_delta[key].append(list[i])

	df_delta = pd.DataFrame.from_dict(dict_delta)

	with pd.ExcelWriter(f"results_gamma_{gamma_folder}.xlsx") as writer:
		df_delta.to_excel(writer, sheet_name="Deltas", index=False)

if gen_curves_flag == True:
	for ax in [ax1,ax2,ax3,ax4]:
		ax.set_xlabel(r'a$_{\text{L}}$',fontsize=22)
		ax.legend(fontsize=16)
	for fig in [fig1,fig2,fig3,fig4]:
		fig.suptitle(f'{to_latex_formula(name_bin)}',fontsize=22, y=0.95)
	ax1.set_ylabel(r'$\Delta$F (k$_{\text{b}}$T)',fontsize=22)
	ax2.set_ylabel(r'$\Delta$U (k$_{\text{b}}$T)',fontsize=22)
	ax3.set_ylabel(r'-T$\Delta$S (k$_{\text{b}}$T)',fontsize=22)
	ax4.set_ylabel(r'$\Delta$F (k$_{\text{b}}$T)',fontsize=22)
	fig1.savefig(f"{final_output}/F_binary.png", format="png", dpi=300,bbox_inches='tight')
	#fig4.savefig(f"{final_output}/F_binary_pairwise.png", format="png", dpi=300,bbox_inches='tight')
	fig2.savefig(f"{final_output}/U_binary.png", format="png", dpi=300,bbox_inches='tight')
	fig3.savefig(f"{final_output}/S_binary.png", format="png", dpi=300,bbox_inches='tight')
	for fig in [fig1,fig2,fig3,fig4]:
		plt.close(fig)

####################### PLOT ###########################
os.chdir(dir_origin)

gamma_value = []
alpha_value = []
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
			gamma = df.loc[df["aL_min"] == "Global X", ["cell"]].iloc[0]
			R2 = df.loc[df["#part"] == "R2 [nm]", ["cell"]].iloc[0]
			alpha = R2/R1_np
			if name != "binary":
				DU, DS, DF = values.sort_values(by="ΔF_min").iloc[0]
			else:
				DU, DS, DF = values.iloc[0]

		elif j == 3:
			values = df.loc[df["aL_min"] == "Global X", ["cell", "ΔU_min", "-ΔST_min", "ΔF_min"]]
			gamma, DU, DS, DF = values.iloc[0]
		
		DU_values[j].append(DU)
		DS_values[j].append(DS)
		DF_values[j].append(DF)

	gamma_value.append(gamma)
	alpha_value.append(alpha)

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
y_label = [r'$\Delta$U (k$_{\text{B}}$T)',r'$-T\Delta$S (k$_{\text{B}}$T)',r'$\Delta$F (k$_{\text{B}}$T)']

ref_pair_data = pd.read_excel(os.path.join(dir_script,"references","data_pair.xlsx"), engine="openpyxl",skiprows=1,header=[0,1])
ref_pair = ref_pair_data[name_bin]
filtered_ref_pair = ref_pair[
    (ref_pair['Gamma'] >= np.min(gamma_value)) & (ref_pair['Gamma'] <= np.max(gamma_value))
]

F_plot = [DU_values[3],DS_values[3],DF_values[3]]
for i, (lista, F) in enumerate(zip(F_plot,["ΔU", "-TΔS", "ΔF"])):
	plt.figure(figsize=(8, 6))
	plt.plot(gamma_value,F_plot[i],ls='none',marker='s',color='red',ms=7,label='MOLT-CF')
	plt.scatter(gamma_MD,F_MD[F],marker='o',color='purple',s=50,label='MD (OTM)',zorder=10)
	plt.scatter(filtered_ref_pair['Gamma'],filtered_ref_pair['DF'],marker='v',color='orange',s=50,label='Pairwise',zorder=10)

	plt.axhline(0,ls='--',c='darkgray',zorder=-3)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.xlim(min(gamma_list)-0.05,max(gamma_list)+0.05)
	plt.ylabel(y_label[i],fontsize=22)
	plt.xlabel(r'$\gamma$',fontsize=22)
	plt.title(f"{to_latex_formula(name_bin)}",fontsize=22,weight='bold')
	plt.legend(fontsize=16)
	plt.savefig(f'{final_output}/{name_bin}_results_{F}.png',format='png',dpi=300,bbox_inches='tight')

F_plot = [DU_values[3],DS_values[3],DF_values[3]]
for i, (lista, F) in enumerate(zip(F_plot,["ΔU", "-TΔS", "ΔF"])):
	plt.figure(figsize=(8, 6))
	plt.plot(alpha_value,F_plot[i],ls='none',marker='s',color='red',ms=7,label='MoltCF')

	plt.axhline(0,ls='--',c='darkgray',zorder=-3)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.xlim(np.min(alpha_value)-0.05,np.max(alpha_value)+0.05)
	plt.ylabel(y_label[i],fontsize=22)
	plt.xlabel(r'$\alpha$',fontsize=22)
	plt.title(f"{to_latex_formula(name_bin)}",fontsize=22,weight='bold')
	plt.legend(fontsize=16)
	plt.savefig(f'{final_output}/{name_bin}_results_{F}_alpha.png',format='png',dpi=300,bbox_inches='tight')


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
		plt.xticks(fontsize=18)
		plt.yticks(fontsize=18)
		plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
		plt.ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=22)
		plt.xlabel(r'$\gamma$',fontsize=22)
		plt.title(f"{to_latex_formula(name_bin)}",fontsize=22,weight='bold')
		plt.legend(fontsize=16)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔU.png',format='png',dpi=300,bbox_inches='tight')

	DS = pd.DataFrame(DS_values).to_numpy()
	plt.figure(figsize=(8, 6))
	for i in range(len(label)):
		plt.plot(gamma_value,DS[i],marker=marker[i],ms=7, label=label[i])

		plt.axhline(0,ls='--',c='darkgray',zorder=-3)
		plt.xticks(fontsize=18)
		plt.yticks(fontsize=18)
		plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
		plt.ylabel(r'$-T\Delta$S (k$_{\text{B}}$T)',fontsize=22)
		plt.xlabel(r'$\gamma$',fontsize=22)
		plt.title(f"{to_latex_formula(name_bin)}",fontsize=22,weight='bold')
		plt.legend(fontsize=16)
		plt.savefig(f'{final_output}/{name_bin}_contrib_results_-TΔS.png',format='png',dpi=300,bbox_inches='tight')

DF = pd.DataFrame(DF_values).to_numpy()
plt.figure(figsize=(8, 6))
for i in range(len(label)):
	plt.plot(gamma_value,DF[i],marker=marker[i],ms=7, label=label[i])

	plt.axhline(0,ls='--',c='darkgray',zorder=-3)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.xlim(min(gamma_value)-0.05,max(gamma_value)+0.05)
	plt.ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=22)
	plt.xlabel(r'$\gamma$',fontsize=22)
	plt.title(f"{to_latex_formula(name_bin)}",fontsize=22,weight='bold')
	plt.legend(fontsize=16)
	plt.savefig(f'{final_output}/{name_bin}_contrib_results_ΔF.png',format='png',dpi=300,bbox_inches='tight')