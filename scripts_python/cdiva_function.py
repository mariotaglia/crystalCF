import subprocess
import csv
import os
import shutil
import pandas as pd
import glob
import re
import numpy as np

def read_DEF(file_path):
    """Extract the lines from DEF."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines

def write_DEF(file_path, lines):
    """Export DEF with lines."""
    with open(file_path, 'w') as f:
        f.writelines(lines)
        
def update_cdiva(DEF, name_bin, gamma, flag_reflexion):
    lines = read_DEF(DEF)
    if name_bin == 'MgZn2':
        for i, line in enumerate(lines):
            if line == "!cdiva\n":
                size_index = i + 1
                cdiva = float(lines[size_index].split()[0])
                lines[size_index] = f"{str(cdiva/2)}\n"
                break

    names = ["CaCu5", "AlB2", "AuCu"]
    if name_bin in names:
        for i, line in enumerate(lines):
            if line == "!cdiva\n":
                size_index = i + 1
                cdiva = cdiva_calc(name_bin,gamma)
                lines[size_index] = f"{str(cdiva)}\n"
                break
    if name_bin == 'NaZn13':
        cdiva, vector = pos_func_gamma(gamma, name_bin)

        for i, line in enumerate(lines):
            if line == "!cdiva\n":
                size_index = i + 1
                lines[size_index] = f"{str(cdiva)}\n"

            if line.strip() == "!Center":
                for j, vec in enumerate(vector):
                    size_index = i + 1 + j
                    pos_list = vec
                    pos1, pos2, pos3 = [pos_list[0],pos_list[1],pos_list[2]]
                    lines[size_index] = f"{pos1} {pos2} {pos3}\n"

    if name_bin == 'CaB6':
        cdiva, vector = pos_func_gamma(gamma, name_bin)

        for i, line in enumerate(lines):
            if line == "!cdiva\n":
                size_index = i + 1
                lines[size_index] = f"{str(cdiva)}\n"

            if line.strip() == "!Center":
                for j, vec in enumerate(vector):
                    size_index = i + 1 + j
                    pos_list = vec
                    pos1, pos2, pos3 = [pos_list[0],pos_list[1],pos_list[2]]
                    lines[size_index] = f"{pos1} {pos2} {pos3}\n"

    write_DEF("DEFINITIONS.txt", lines)

def cdiva_calc(name,gamma):
    if name == "CaCu5":
        if gamma < 2.0 / np.sqrt(3.0) - 1:
            a_fac = 1.0
            c_fac = 1.0
        elif (gamma >= 2.0 / np.sqrt(3.0) - 1) and (gamma < (1 + 2 * np.sqrt(19)) / 15.0):
            a_fac = np.sqrt(3) * (1 + gamma) * 0.5
            c_fac = 1.0
        elif (gamma >= (1 + 2 * np.sqrt(19)) / 15.0) and (gamma<1 / (4.0/np.sqrt(3)-1.0)):
            a_fac = np.sqrt(3) * (1 + gamma)*0.5
            c_fac = 0.5 * np.sqrt(15.0*gamma**2 - 2*gamma - 1)
        else:
            a_fac = 2*gamma
            c_fac = np.sqrt(8.0/3.0)*gamma

    if name == "AlB2":
        if gamma < np.sqrt(7.0/3.0)-1:
            a_fac = 1.0
            c_fac = 1.0
        elif gamma > 2.0/3.0:
            a_fac = np.sqrt(3.0)*gamma
            c_fac = 1.0
        elif (gamma > 1/np.sqrt(3.0)) and (gamma <= 2.0/3.0):
            a_fac = np.sqrt(3.0)*gamma
            c_fac = np.sqrt(-3.0*gamma**2+2.0*gamma+1.0)
        else:
            c_fac = np.sqrt(gamma**2+2*gamma-1.0/3.0)
            a_fac = 1.0

    if name == "AuCu":
        if (0 < gamma) and (gamma < np.sqrt(3.0)-1):
            a_fac = np.sqrt(2)
            c_fac = np.sqrt(2)

        elif (gamma >= np.sqrt(3.0)-1) and (gamma <= 1.0):
            a_fac = np.sqrt(2)
            c_fac = np.sqrt(2)*np.sqrt(0.5*gamma**2+gamma-0.5)

        else:
            raise ValueError('Parameter gamma is outside allowed range')

    if name == "Cu3Au":
        a_fac = 1
        if gamma < np.sqrt(2.0)-1:
            c_fac = 1
        else:
            c_fac = (1+gamma)/np.sqrt(2.0)

    cdiva = c_fac/a_fac
    if cdiva >= 2 or cdiva <= 0.5:
        print(f"warning, cdiva<=0.5 or cdiva>=2.0 for gamma {gamma:.2f}")
    return cdiva

def pos_func_gamma(gam, name_bin):
    if name_bin == "NaZn13":
        a_nn_e = 1.0  # distancia efectiva entre primeros vecinos
        l_value = 1.0  # escala general (usualmente 1.0)

        # Cálculos preliminares
        tau = 0.5 * (1 + np.sqrt(5.0))
        theta = 1 + 2 * (1 + tau) / (np.sqrt(1 + tau ** 2))
        psi = 2 * (11.0 + np.sqrt(5) + 2 * np.sqrt(10.0 * (1 + np.sqrt(5)))) / (5 + np.sqrt(5))
        gamma_c1 = (theta - np.sqrt(theta ** 2 - 6.0)) / 3.0
        gamma_c2 = (1 + np.sqrt(psi + 1)) / psi

        # Factores de escala
        if gam < gamma_c1:
          c_fac = 2.0
        elif (gam >= gamma_c1) and (gam < gamma_c2):
          t_fac = (3.0 + np.sqrt(5)) * np.sqrt(2.0 / (5.0 + np.sqrt(5.0)))
          c_fac = 2.0 * (t_fac * gam + np.sqrt(3 + 6 * gam + (8.0 / np.sqrt(5) - 5) * gam ** 2)) / 3.0
        else:
          c_fac = 2 * gam * (2 * tau + np.sqrt(tau ** 2 - 1)) / np.sqrt(1 + tau ** 2)

        cdiva = 1.0
        # Vectores de red
        a_axis = b_axis = c_axis = l_value * c_fac
        a_vector = np.zeros((3, 3))
        a_vector[0] = a_axis * np.array([1.0, 0.0, 0.0])
        a_vector[1] = b_axis * np.array([0.0, 1.0, 0.0])
        a_vector[2] = c_axis * np.array([0.0, 0.0, 1.0])

        # Posiciones de tipo A (8 átomos)
        hb = np.zeros((8, 3))
        hb[1] = [0.0, 0.5, 0.5]
        hb[2] = [0.5, 0.5, 0.0]
        hb[3] = [0.5, 0.0, 0.5]
        hb[4] = [0.5, 0.5, 0.5]
        hb[5] = [0.5, 0.0, 0.0]
        hb[6] = [0.0, 0.5, 0.0]
        hb[7] = [0.0, 0.0, 0.5]

        v_vector = []

        # Tipo A
        v_vec = np.array([0.25, 0.25, 0.25])
        for i in range(8):
          v_vector.append((v_vec + hb[i]) @ a_vector)

        # Primeros icosaedros
        modt = np.sqrt(1 + tau ** 2)
        ico = np.zeros((13, 3))
        ico[1] = [0, 1, tau]
        ico[2] = [0, 1, -tau]
        ico[3] = [0, -1, tau]
        ico[4] = [0, -1, -tau]
        for i in range(4):
          ico[5 + i] = np.roll(ico[i + 1], 1)
          ico[9 + i] = np.roll(ico[i + 1], 2)

        for i in range(4):
          for j in range(13):
              pos = a_nn_e * gam * ico[j] / modt + hb[i] @ a_vector
              v_vector.append(pos)

        # Segundos icosaedros
        ico = np.zeros((13, 3))
        ico[1] = [1, 0, tau]
        ico[2] = [-1, 0, tau]
        ico[3] = [1, 0, -tau]
        ico[4] = [-1, 0, -tau]
        for i in range(4):
          ico[5 + i] = np.roll(ico[i + 1], 1)
          ico[9 + i] = np.roll(ico[i + 1], 2)

        for i in range(4):
          for j in range(13):
              pos = a_nn_e * gam * ico[j] / modt + hb[4 + i] @ a_vector
              v_vector.append(pos)

        # Convertir a array y normalizar en fracción de celda
        v_vector = np.array(v_vector)
        vector_primer_oct = v_vector[np.all((v_vector/c_fac >= 0.0) & (v_vector/c_fac <= 0.5), axis=1)]*2/c_fac
        return cdiva, np.round(vector_primer_oct,10)

    if name_bin == "CaB6":
        if gam < 1 / (1 + np.sqrt(2)):
            c_fac = 1
            u_val = (1 - np.sqrt(2) * gam) / 2.0
        else:
            c_fac = gam * (1 + np.sqrt(2))
            u_val = 1 / (2.0 * np.sqrt(2) + 2.0)

        cdiva = 1.0
        v_vector = []

        v_vector.append(np.array([0.0, 0.0, 0.0]))
        v_vector.append(np.array([0.5, 0.5, u_val]))
        v_vector.append(np.array([0.5, 0.5, 1.0 - u_val]))
        v_vector.append(np.array([u_val, 0.5, 0.5]))
        v_vector.append(np.array([1 - u_val, 0.5, 0.5]))
        v_vector.append(np.array([0.5, u_val, 0.5]))
        v_vector.append(np.array([0.5, 1 - u_val, 0.5]))
        v_vector = np.array(v_vector)

        return cdiva, np.round(v_vector,10)