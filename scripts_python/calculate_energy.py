import os
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import math
from ast import literal_eval

def mean_al(df):
    return df.groupby('aL', as_index=False).agg({
        'aL': 'mean',
        'F_tot_gcanon': 'mean',
        'F_norm': 'mean'
    })

def mean_al_eta(df):
    return df.groupby('aL', as_index=False).agg({
        'aL': 'mean',
        'eta': 'mean',
    })

def mean_al_vol(df):
    return df.groupby("aL")["Vol"].apply(lambda x: np.mean(np.stack(x), axis=0).tolist()).reset_index()


def estimate_part_F(part, part_cell, factor_aL_part, ni, k_part, n1, n2, gen_curves_flag, k_reflex_part):
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

    F_min_cell = y_cell.min()/(n1+n2)
    aL_min_cell = x_cell[y_cell.argmin()]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_norm_cell/(n1+n2))
        ax_sub.plot(x_cell,y_cell/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_cell,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_{part}_{part_cell}.png", format="png")
        plt.close(fig_sub)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_tot/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_cell,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_without_ref_{part}_{part_cell}.png", format="png")
        plt.close(fig_sub)

    return aL_min_cell, F_min_cell, x_cell

def estimate_part_contrib(part, part_cell, factor_aL_part, ni, k_part, n1, n2, F, aL_array, aL_min, gen_curves_flag, k_reflex_part):
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
    F_min_cell = y_cell[aL_array==aL_min]/(n1+n2)
    F_calc = F_min_cell[0]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_cell, F_norm_cell/(n1+n2))
        ax_sub.plot(aL_array,y_cell/(n1+n2))
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        if "F_trans" in F:
            ax_sub.set_ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"U_{part}_{part_cell}.png", format="png")
        if "F_HS" in F:
            ax_sub.set_ylabel(r'-T$\Delta$S (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"S_{part}_{part_cell}.png", format="png")
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

    F_min_bin = y_bin.min()/(n1+n2)
    aL_min_bin = x_bin[y_bin.argmin()]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_norm_bin/(n1+n2))
        ax_sub.plot(x_bin,y_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_bin,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_{name}.png", format="png")
        plt.close(fig_sub)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_tot_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\Delta$F (k$_{\text{B}}$T)',fontsize=14)
        ax_sub.axvline(aL_min_bin,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"F_without_ref_{name}.png", format="png")
        plt.close(fig_sub)

        ax.scatter(aL_bin, F_norm_bin/(n1+n2), label=fr'$\gamma:$ {gamma}')
        ax.plot(x_bin,y_bin/(n1+n2))
        ax.scatter(aL_min_bin,F_min_bin,marker='|', color='black',s=50,zorder=10)

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

    F_min_bin = y_bin[aL_array==aL_min]/(n1+n2)
    F_calc = F_min_bin[0]

    if gen_curves_flag == True:
        ax.scatter(aL_bin, F_norm_bin/(n1+n2), label=fr'$\gamma:$ {gamma}')
        ax.plot(aL_array,y_bin/(n1+n2))
        ax.scatter(aL_min,F_min_bin,marker='|', color='black',s=50,zorder=10)

        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, F_norm_bin/(n1+n2))
        ax_sub.plot(aL_array,y_bin/(n1+n2))
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        if "F_trans" in F:
            ax_sub.set_ylabel(r'$\Delta$U (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"U_{name}.png", format="png")
        if "F_HS" in F:
            ax_sub.set_ylabel(r'-T$\Delta$S (k$_{\text{B}}$T)',fontsize=14)
            fig_sub.savefig(f"S_{name}.png", format="png")
        plt.close(fig_sub)

    return F_calc

def delta_energy_F(dict_delta, cell, n1, n2):
    df_delta = pd.DataFrame.from_dict(dict_delta)
    F_part1 = (df_delta.loc[df_delta["#part"] == "part1", "ΔF_min"]*n1).tolist()
    F_part2 = (df_delta.loc[df_delta["#part"] == "part2", "ΔF_min"]*n2).tolist()
    F_min_bin = df_delta.loc[df_delta["#part"] == "binary", "ΔF_min"].values[0]  # Tomar único valor binario
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
    F_part1 = df_delta.loc[(df_delta["#part"] == "part1") & (df_delta["cell"] == cell_min[0]), f"{US}_min"].values[0]
    F_part2 = df_delta.loc[(df_delta["#part"] == "part2") & (df_delta["cell"] == cell_min[1]), f"{US}_min"].values[0]
    F_min_bin = df_delta.loc[df_delta["#part"] == "binary", f"{US}_min"].values[0]  # Tomar único valor binario
    DF = F_min_bin - F_part1*n1 - F_part2*n2

    return DF

def delta_energy_contrib(dict_contrib, File_name, cell_min, n1, n2):
    DF_calc = []
    df_contrib = pd.DataFrame.from_dict(dict_contrib)
    for i, F in enumerate(File_name):
        F_part1 = df_contrib.loc[(df_contrib["#part"] == "part1") & (df_contrib["cell"] == cell_min[0]), f"{F}"].values[0]
        F_part2 = df_contrib.loc[(df_contrib["#part"] == "part2") & (df_contrib["cell"] == cell_min[1]), f"{F}"].values[0]
        F_min_bin = df_contrib.loc[df_contrib["#part"] == "binary", f"{F}"].values[0]  # Tomar único valor binario

        DF_calc.append(F_min_bin - F_part1*n1 - F_part2*n2)

    return DF_calc

def estimate_part_volume_overlap(part, part_cell, factor_aL_part, k_part, aL_array, aL_min, k_reflex_part, gen_curves_flag, ax1, ax2):
    csv_file = f"{part}_results_volume_overlap.csv"
    data_part = pd.read_csv(csv_file, skiprows=0, delimiter=';')

    data_part_cell = data_part[data_part["cell"] == part_cell].copy()
    data_part_cell['aL'] = data_part_cell['delta'] * data_part_cell['dimx'] * factor_aL_part*k_reflex_part
    data_part_cell['aL'] = data_part_cell['aL'].round(4)
    data_part_cell["Vol"] = data_part_cell['volume overlap [nm³]'].apply(lambda x: [float(i) for i in literal_eval(x)])
    data_part_cell.sort_values(by='aL', inplace=True)

    k_reflex = k_reflex_part**3
    df_cell = mean_al_vol(data_part_cell)
    aL_cell = df_cell['aL'].to_numpy()
    
    V_cell = np.array([np.array(x) for x in df_cell['Vol']])
    V_sum = df_cell['Vol'].apply(lambda x: np.sum(np.array(x), axis=0)).to_numpy()/k_part

    y_cell = []; V_cell_min = []
    for i in range(2):
        cs = CubicSpline(aL_cell, V_cell.T[i])(aL_array)
        y_cell.append(cs)
        V_cell_min.append(cs[aL_array==aL_min])

    y_sum_cell = CubicSpline(aL_cell, V_sum)(aL_array)
    V_sum_cell_min = y_sum_cell[aL_array==aL_min]

    if gen_curves_flag == True:
        for i in range(2):
            ax1.scatter(aL_cell, V_cell.T[i], label=part_cell)
            ax1.plot(aL_array,y_cell[i])
            ax1.scatter(aL_min,V_cell_min[i],marker='|', color='black',s=50,zorder=10)

        ax2.scatter(aL_cell, V_sum, label=part_cell)
        ax2.plot(aL_array,y_sum_cell)
        ax2.scatter(aL_min,V_sum_cell_min,marker='|', color='black',s=50,zorder=10)
        for ax in [ax1,ax2]:
            ax.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
            ax.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
            ax.set_ylabel(r'Vol. overlap [nm$^3$]',fontsize=14)
            ax.legend()

    return V_cell_min, V_sum_cell_min

def estimate_bin_volume_overlap(name, factor_bin_cell, k_bin, k_reflex, n1, n2, aL_array, aL_min, gamma, gen_curves_flag):
    V_bin_min = {"part1": [], "part2": []}
    k = k_reflex["kx"]*k_reflex["ky"]*k_reflex["kz"]

    csv_file = f"{name}_results_volume_overlap.csv"
    data_bin = pd.read_csv(csv_file, skiprows=0, delimiter=';')

    data_bin['aL'] = (data_bin['delta'] * data_bin['dimx'] * float(factor_bin_cell)*k_reflex["kx"]).round(4)
    data_bin["Vol"] = data_bin['volume overlap [nm³]'].apply(lambda x: [float(i) for i in literal_eval(x)])
    
    # Ordenar por 'aL'
    data_bin.sort_values(by='aL', inplace=True)
    aL_bin = data_bin['aL'].to_numpy()

    df_bin = mean_al_vol(data_bin)
    aL_bin = df_bin['aL'].to_numpy()
    
    V_bin = np.array([np.array(x) for x in df_bin['Vol']])
    V_sum = df_bin['Vol'].apply(lambda x: np.sum(np.array(x), axis=0)).to_numpy()/k_bin

    y_bin = []; V_min = []
    for i in range((n1+n2)*k_bin):
        cs = CubicSpline(aL_bin, V_bin.T[i])(aL_array)
        y_bin.append(cs)
        V_min.append(cs[aL_array==aL_min])

    for j, part in zip([(0,n1*k_bin),(n1*k_bin,(n1+n2)*k_bin)],["part1", "part2"]):
        for k in range(j[0],j[1]):
            V_bin_min[part].append(V_min[k])

    y_sum_bin = CubicSpline(aL_bin, V_sum)(aL_array)
    V_sum_bin_min = y_sum_bin[aL_array==aL_min]

    if gen_curves_flag == True:
        fig, (ax1,ax2) = plt.subplots(ncols = 2, nrows = 1)
        for i in range((n1+n2)*k_bin):
            ax1.scatter(aL_bin, V_bin.T[i], label=fr'$part {i}')
            ax1.plot(aL_array,y_bin[i])
            ax1.scatter(aL_min,V_min[i],marker='|', color='black',s=50,zorder=10)

        ax2.scatter(aL_bin, V_sum, label=fr'$\gamma:$ {gamma}')
        ax2.plot(aL_array,y_sum_bin)
        ax2.scatter(aL_min,V_sum_bin_min,marker='|', color='black',s=50,zorder=10)
        for ax in [ax1,ax2]:
            ax.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
            ax.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
            ax.set_ylabel(r'Vol. overlap [nm$^3$]',fontsize=14)
            ax.legend()

        fig.savefig(f"Vol_overlap_{name}.png", format="png")
        plt.close(fig)

    return V_bin_min, V_sum_bin_min

def packing_frac_part(part, part_cell, factor_aL_part, k_part, aL_array, aL_min, k_reflex_part, gen_curves_flag):
    csv_file = f"{part}_results_volume_overlap.csv"
    data_part = pd.read_csv(csv_file, skiprows=0, delimiter=';')

    data_part_cell = data_part[data_part["cell"] == part_cell].copy()
    data_part_cell['aL'] = data_part_cell['delta'] * data_part_cell['dimx'] * factor_aL_part*k_reflex_part
    data_part_cell['aL'] = data_part_cell['aL'].round(4)
    data_part_cell["eta"] = data_part_cell['packing fraction']
    data_part_cell.sort_values(by='aL', inplace=True)

    df_part = mean_al_eta(data_part_cell)
    aL_part = df_part['aL'].to_numpy()
    eta_part = df_part["eta"].to_numpy()

    y_cell = CubicSpline(aL_part, eta_part)(aL_array)
    eta_min = y_cell[aL_array==aL_min]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_part, eta_part)
        ax_sub.plot(aL_array,y_cell)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\eta$',fontsize=14)
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"eta_{part}_{part_cell}.png", format="png")
        plt.close(fig_sub)
    return eta_min[0]

def packing_frac_bin(name, factor_bin_cell, aL_array, aL_min, k_reflex, gen_curves_flag):
    csv_file = f"{name}_results_volume_overlap.csv"
    data_bin = pd.read_csv(csv_file, skiprows=0, delimiter=';')
    k = k_reflex["kx"]*k_reflex["ky"]*k_reflex["kz"]
    data_bin['aL'] = (data_bin['delta'] * data_bin['dimx'] * float(factor_bin_cell)*k_reflex["kx"]).round(4)
    data_bin["eta"] = data_bin['packing fraction']
    data_bin.sort_values(by='aL', inplace=True)

    df_bin = mean_al_eta(data_bin)
    aL_bin = df_bin['aL'].to_numpy()
    eta_bin = df_bin["eta"].to_numpy()

    y_bin = CubicSpline(aL_bin, eta_bin)(aL_array)
    eta_min = y_bin[aL_array==aL_min]

    if gen_curves_flag == True:
        fig_sub, ax_sub = plt.subplots()
        ax_sub.scatter(aL_bin, eta_bin)
        ax_sub.plot(aL_array,y_bin)
        ax_sub.set_xlabel(r'a$_{\text{L}}$',fontsize=14)
        ax_sub.set_ylabel(r'$\eta$',fontsize=14)
        ax_sub.axvline(aL_min,ls='--',c='darkgray',zorder=-1)
        fig_sub.savefig(f"eta_{name}.png", format="png")
        plt.close(fig_sub)
    return eta_min[0]
