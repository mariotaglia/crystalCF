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

    names = ["CaCu5", "AlB2"]
    if name_bin in names:
        for i, line in enumerate(lines):
            if line == "!cdiva\n":
                size_index = i + 1
                cdiva = cdiva_calc(name_bin,gamma)
                lines[size_index] = f"{str(cdiva)}\n"
                break

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

    cdiva = c_fac/a_fac
    if cdiva > 2 or cdiva <= 0.5:
        print(f"warning, cdiva<0.5 or cdiva>2.0 for gamma {gamma:.2f}")
    return cdiva