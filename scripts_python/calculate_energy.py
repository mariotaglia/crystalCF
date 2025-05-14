import os
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import math

def mean_al(df):
    return df.groupby('aL', as_index=False).agg({
        'aL': 'mean',
        'F_tot_gcanon': 'mean',
        'F_norm': 'mean'
    })

def estimate_part_F(part, part_cell, factor_aL_part, ni, k_part, gen_curves_flag, k_reflex_part):
    F_norm = []
    csv_file = [f"{part}_results_output.csv", f"{part}_references_output.csv"]
    data_part = pd.read_csv(csv_file[0], skiprows=0)
    data_part_ref = pd.read_csv(csv_file[1], skiprows=0)
    data_part_cell = data_part[data_part["cell"] == part_cell].copy()
    data_part_cell['aL'] = data_part_cell['delta'] * data_part_cell['dimx'] *factor_aL_part*k_reflex_part
    data_part_cell['aL'] = data_part_cell['aL'].round(4) #needed to calculate mean.
    data_part_cell["F_norm"] = data_part_cell["F_tot_gcanon"] - data_part_ref["F_tot_gcanon_reference"]
    data_part_cell.sort_values(by='aL', inplace=True)
    df_cell = mean_al(data_part_cell)
    df_tot_cell = mean_al(data_part_cell)

    k_reflex = k_reflex_part**3
    aL_cell = df_cell['aL'].to_numpy()
    F_tot = df_cell["F_tot_gcanon"].to_numpy()*k_reflex/k_part
    F_norm_cell = df_cell['F_norm'].to_numpy()*k_reflex/k_part

    x_cell =  np.arange(aL_cell[0], aL_cell[-1], 0.001)
    y_cell = CubicSpline(aL_cell, F_norm_cell)(x_cell)

    # Ajuste cúbico (polinomio de grado 3)
    coeficientes = np.polyfit(aL_cell, F_norm_cell, 4)
    polinomio = np.poly1d(coeficientes)

    # Evaluación del polinomio ajustado
    y_cell = polinomio(x_cell)
    F_min_cell = y_cell.min()
    aL_min_cell = x_cell[y_cell.argmin()]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_norm_cell)
        ax_sub.plot(x_cell,y_cell)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_cell,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_{part}_{part_cell}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_tot)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_cell,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_without_ref_{part}_{part_cell}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

    return aL_min_cell, F_min_cell, x_cell

def estimate_part_contrib(part, part_cell, factor_aL_part, ni, k_part, F, aL_array, aL_min, gen_curves_flag, k_reflex_part):
    F_cell = []
    for f_name in F:
        csv_file = [f"{part}_results_output.csv", f"{part}_references_output.csv"]

        data_part = pd.read_csv(csv_file [0], skiprows=0)
        data_part_ref = pd.read_csv(csv_file [1], skiprows=0)

        data_part_cell = data_part[data_part["cell"] == part_cell].copy()
        data_part_cell['aL'] = data_part_cell['delta'] * data_part_cell['dimx'] * factor_aL_part*k_reflex_part
        data_part_cell['aL'] = data_part_cell['aL'].round(4) #si no lo redondeo no puede promediar.
        data_part_cell["F_norm"] = data_part_cell[f_name] - data_part_ref[f"{f_name}_reference"]
        data_part_cell.sort_values(by='aL', inplace=True)

        k_reflex = k_reflex_part**3
        df_cell = mean_al(data_part_cell)
        aL_cell = df_cell['aL'].to_numpy()
        F_cell.append(df_cell['F_norm'].to_numpy()*k_reflex/k_part)

    F_norm_cell = np.sum(F_cell, axis=0)
    y_cell = CubicSpline(aL_cell, F_norm_cell)(aL_array)
    
    coeficientes = np.polyfit(aL_cell, F_norm_cell, 4)
    polinomio = np.poly1d(coeficientes)
    y_cell = polinomio(aL_array)

    F_min_cell = y_cell[aL_array==aL_min]
    F_calc = F_min_cell[0]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_norm_cell)
        ax_sub.plot(aL_array,y_cell)
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        if "F_trans" in F:
            ax_sub.set_ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"U_{part}_{part_cell}.png", format="png",bbox_inches='tight')
        if "F_HS" in F:
            ax_sub.set_ylabel(r'-T$\Delta$S (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"S_{part}_{part_cell}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

    return F_calc

def estimate_bin_F(name, factor_bcell, k_bin, n1, n2, ax, gamma, gen_curves_flag, k_reflex):
    F_norm = []
    csv_file = [f"{name}_results_output.csv", f"{name}_references_output.csv"]
    data_bin = pd.read_csv(csv_file[0], skiprows=0)
    data_bin_ref = pd.read_csv(csv_file[1], skiprows=0)

    k = k_reflex["kx"]*k_reflex["ky"]*k_reflex["kz"]
    for part in ["part1", "part2"]:
        data_bin_part = data_bin_ref[data_bin_ref["#part"] == part].copy()
        data_bin_part = data_bin_part.rename(columns={"F_tot_gcanon_reference": f"F_tot_gcanon_reference_{part}"})
        data_bin = data_bin.merge(data_bin_part[['delta', 'dimx', f'F_tot_gcanon_reference_{part}']], 
                                  on=['delta', 'dimx'], 
                                  how='left')
        
    data_bin['aL'] = (data_bin['delta'] * data_bin['dimx'] * float(factor_bcell)*k_reflex["kx"]).round(4)
    data_bin["F_norm"] = data_bin["F_tot_gcanon"] - data_bin["F_tot_gcanon_reference_part1"] - data_bin["F_tot_gcanon_reference_part2"]
    data_bin.sort_values(by='aL', inplace=True)

    df_bin = mean_al(data_bin)
    aL_bin = df_bin['aL'].to_numpy()
    F_norm_bin = df_bin['F_norm'].to_numpy()*k/k_bin
    F_tot_bin = df_bin["F_tot_gcanon"].to_numpy()*k/k_bin
    x_bin =  np.arange(aL_bin[0], aL_bin[-1], 0.001)
    y_bin = CubicSpline(aL_bin, F_norm_bin)(x_bin)

    coeficientes = np.polyfit(aL_bin, F_norm_bin, 4)
    polinomio = np.poly1d(coeficientes)
    y_bin = polinomio(x_bin)

    F_min_bin = y_bin.min()
    aL_min_bin = x_bin[y_bin.argmin()]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_norm_bin/(n1+n2))
        ax_sub.plot(x_bin,y_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_bin,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_{name}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_tot_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_bin,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_without_ref_{name}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

        ax.scatter(aL_bin, F_norm_bin/(n1+n2), label=fr'$\gamma:$ {gamma}')
        ax.plot(x_bin,y_bin/(n1+n2))
        ax.scatter(aL_min_bin,F_min_bin/(n1+n2),marker='|', color='black',s=50,zorder=10)

    return aL_min_bin, F_min_bin, x_bin

def estimate_bin_contrib(name, factor_bin_cell, k_bin, n1, n2, F, aL_array, aL_min, ax, gamma, gen_curves_flag, k_reflex):
    F_bin = []
    k = k_reflex["kx"]*k_reflex["ky"]*k_reflex["kz"]
    for f_name in F:
        csv_file = [f"{name}_results_output.csv", f"{name}_references_output.csv"]
        data_bin = pd.read_csv(csv_file[0], skiprows=0)
        data_bin_ref = pd.read_csv(csv_file[1], skiprows=0)

        for part in ["part1", "part2"]:
            data_bin_part = data_bin_ref[data_bin_ref["#part"] == part].copy()
            data_bin_part = data_bin_part.rename(columns={f"{f_name}_reference": f"{f_name}_reference_{part}"})
            data_bin = data_bin.merge(data_bin_part[['delta', 'dimx', f"{f_name}_reference_{part}"]], 
                                      on=['delta', 'dimx'], 
                                      how='left')

        data_bin['aL'] = (data_bin['delta'] * data_bin['dimx'] * float(factor_bin_cell)*k_reflex["kx"]).round(4)
        data_bin["F_norm"] = data_bin[f_name] - data_bin[f"{f_name}_reference_part1"] - data_bin[f"{f_name}_reference_part2"]

        # Ordenar por 'aL'
        data_bin.sort_values(by='aL', inplace=True)

        # Promediar valores repetidos de aL
        df_bin = mean_al(data_bin)
        aL_bin = df_bin['aL'].to_numpy()
        F_bin.append(df_bin['F_norm'].to_numpy()*k/k_bin)

    F_norm_bin = np.sum(F_bin, axis=0)
    y_bin = CubicSpline(aL_bin, F_norm_bin)(aL_array)

    coeficientes = np.polyfit(aL_bin, F_norm_bin, 4)
    polinomio = np.poly1d(coeficientes)
    y_bin = polinomio(aL_array)

    F_min_bin = y_bin[aL_array==aL_min]
    F_calc = F_min_bin[0]

    if gen_curves_flag == True:
        ax.scatter(aL_bin, F_norm_bin/(n1+n2), label=fr'$\gamma:$ {gamma}')
        ax.plot(aL_array,y_bin/(n1+n2))
        ax.scatter(aL_min,F_min_bin/(n1+n2),marker='|', color='black',s=50,zorder=10)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_norm_bin/(n1+n2))
        ax_sub.plot(aL_array,y_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        if "F_trans" in F:
            ax_sub.set_ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"U_{name}.png", format="png",bbox_inches='tight')
        if "F_HS" in F:
            ax_sub.set_ylabel(r'-T$\Delta$S (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"S_{name}.png", format="png",bbox_inches='tight')
        plt.close(fig_sub)

    return F_calc

def delta_energy_F(dict_delta, cell, n1, n2):
    df_delta = pd.DataFrame.from_dict(dict_delta)
    F_part1 = (df_delta.loc[df_delta["#part"] == "part1", "ΔF_min"]*n1/(n1+n2)).tolist()
    F_part2 = (df_delta.loc[df_delta["#part"] == "part2", "ΔF_min"]*n2/(n1+n2)).tolist()
    F_min_bin = df_delta.loc[df_delta["#part"] == "binary", "ΔF_min"].values[0]/(n1+n2)  # Tomar único valor binario
    struc_part_min = []
    for i, part in enumerate(["part1", "part2"]):
        F_part = F_part1 if part == "part1" else F_part2  # Seleccionar la lista correcta
        valores = {cell[j]: F_part[j] for j in range(len(F_part))}  # Construir el diccionario correctamente
        
        min_pos = min(valores, key=valores.get)  # Encontrar la estructura con menor ΔF_min
        min_part_value = valores[min_pos]  # Obtener el valor mínimo
        struc_part_min.append(min_pos)
        # Calcular DF correctamente
        if i == 0:  # Primer ciclo
            DF = F_min_bin - min_part_value
            min_part_sum = min_part_value
        else:  # Resto de iteraciones
            DF -= min_part_value
            min_part_sum += min_part_value
    return DF, struc_part_min

def delta_energy_US(dict_delta, US, cell_min, n1, n2):
    df_delta = pd.DataFrame.from_dict(dict_delta)
    F_part1 = df_delta.loc[(df_delta["#part"] == "part1") & (df_delta["cell"] == cell_min[0]), f"{US}_min"].values[0]/(n1+n2)
    F_part2 = df_delta.loc[(df_delta["#part"] == "part2") & (df_delta["cell"] == cell_min[1]), f"{US}_min"].values[0]/(n1+n2)
    F_min_bin = df_delta.loc[df_delta["#part"] == "binary", f"{US}_min"].values[0]/(n1+n2)  # Tomar único valor binario
    DF = F_min_bin - F_part1*n1 - F_part2*n2

    return DF

def delta_energy_contrib(dict_contrib, File_name, cell_min, n1, n2):
    DF_calc = []
    df_contrib = pd.DataFrame.from_dict(dict_contrib)
    for i, F in enumerate(File_name):
        F_part1 = df_contrib.loc[(df_contrib["#part"] == "part1") & (df_contrib["cell"] == cell_min[0]), f"{F}"].values[0]/(n1+n2)
        F_part2 = df_contrib.loc[(df_contrib["#part"] == "part2") & (df_contrib["cell"] == cell_min[1]), f"{F}"].values[0]/(n1+n2)
        F_min_bin = df_contrib.loc[df_contrib["#part"] == "binary", f"{F}"].values[0]/(n1+n2)  # Tomar único valor binario

        DF_calc.append(F_min_bin - F_part1*n1 - F_part2*n2)

    return DF_calc

