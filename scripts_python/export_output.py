import os
import subprocess
import shutil
import numpy as np
import pandas as pd
import csv
from collections import defaultdict

def run_command(command):
    """Run the command on the terminal and return it's output."""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    return result.stdout.strip()

def process_principal(output_file, name_bin, R, delta_dim_bin, aL, k_aL, F):
    structure = os.getcwd()
    R1_np = R["part1"]; R2_np = R["part2"]
    k = 1
    if name_bin == "MgZn2":
        k = 2
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("radius A [nm],radius B [nm],delta,dimx,dimy,dimz,F_value\n")

    delta_list = sorted({entry["delta"] for entry in delta_dim_bin if entry["delta"] is not None})
    for delta in delta_list:
        round_value = int(np.round(float(aL/k_aL["kx"]) / float(delta)))
        dims = []
        dims_sum_bin = [entry["dim"] for entry in delta_dim_bin if entry["delta"] == delta][0]
        #dims_sum_bin = dims_sum_bin = np.arange(-8,16,1)
        for sum_dim in dims_sum_bin:
            dims.append(round_value + int(sum_dim))
        delta_folder = str(delta).replace('.', '_')

        for dim in dims:
            folder_name = f"delta_{delta_folder}_dim_{dim}"
            os.chdir(folder_name)
            
            file_path = os.path.join(os.getcwd(), F+".dat")

            if os.path.isfile(file_path):
                try:
                    with open(file_path, "r") as file:
                        first_line = file.readline().strip().split()
                        if len(first_line) > 1:
                            lista = ['Li3Bi','NaZn13', 'bccAB6']
                            #if name_bin in lista:
                            #    value = float(first_line[1])/8
                            #else:
                            value = first_line[1]

                            with open(output_file, "a") as out_file:
                                out_file.write(f"{R1_np},{R2_np},{delta},{dim},{int(dim*k_aL['kx']/k_aL['ky'])},{int(dim*k_aL['kx']*k/k_aL['kz'])},{value}\n")
                        else:
                            print(f"Advertencia: No se encontró un valor en {file_path}")
                except Exception as e:
                    print(f"Error al procesar {file_path}: {e}")
            else:
                print(f"Advertencia: Archivo no encontrado en {file_path}")
            os.chdir("..")

def process_reference_bin(output_file, dir_inicial, F, R, gamma_folder):
    gamma = gamma_folder.replace('_', '.')
    dir_fuente = {"part1": os.path.join(dir_inicial,"sim_part1"),"part2": os.path.join(os.getcwd())}
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

    for labels in ["part1", "part2"]:
        references_file = os.path.join(dir_fuente[labels], "binary_ref", "tot_references.csv")
        contadores = defaultdict(lambda: defaultdict(int))
        delta_map = defaultdict(lambda: defaultdict(set))

        # Leer el archivo de referencias y construir los mapas
        with open(references_file, "r", encoding="utf-8") as ref_file:
            reader = csv.reader(ref_file, delimiter=',')
            next(reader)  # Saltar encabezado
            for parts in reader:
                label = parts[0].strip()
                contador = int(parts[2].strip())
                pos_tuple = (parts[3], parts[4], parts[5])
                delta = parts[6].strip().replace(",", ".")
                dimx = parts[7].strip().replace(",", ".")
                dimy = parts[8].strip().replace(",", ".")
                dimz = parts[9].strip().replace(",", ".")
                key = parts[10].strip()
                dim = (dimx, dimy, dimz)

                # Formatear los valores de gamma a 3 decimales
                gamma_field = parts[11].strip()
                
                # Aceptar si viene una lista (con comas) o un único valor
                gamma_values = gamma_field.split(',') if ',' in gamma_field else [gamma_field]
                gamma_list = [f"{float(g.strip()):.3f}" for g in gamma_values]

                if gamma not in gamma_list:
                    continue

                contadores[delta][(key, dim)] = contador
                delta_map[(label, delta)][key].add(dim)


        # Crear archivo de salida si no existe
        if not os.path.isfile(output_file):
            with open(output_file, "w") as out_file:
                out_file.write("#part,radius [nm],delta,dimx,dimy,dimz,F_reference\n")

        for (label, delta), key_map in delta_map.items():
            base_folder = os.path.join(dir_fuente[label], "binary_ref", label)
            delta_name = delta.replace('.', '_')
            delta_folder = os.path.join(base_folder, f"delta_{delta_name}")

            for key, dims in key_map.items():
                key_folder = os.path.join(delta_folder, key)
                
                # Verificar si el archivo `F_tot_gcanon.dat` existe
                file_path = os.path.join(key_folder, F + ".dat")
                if os.path.isfile(file_path):
                    try:
                        with open(file_path, "r") as file:
                            value = float(file.readline().strip().split()[1])  # Segunda columna
                            for dim in set(dims):
                                multiplicador = contadores.get(delta, {}).get((key, dim), 1)
                                data[label][delta][dim] += value * multiplicador

                    except Exception as e:
                        print(f"Error al procesar {file_path}: {e}")
                else:
                    print(f"Advertencia: Archivo no encontrado en {file_path}")

    for label in ["part1", "part2"]:
        if label in data:  # Verificar si hay datos para esa etiqueta
            with open(output_file, "a", newline="") as csvfile:
                writer = csv.writer(csvfile)
                for delta, dim_values in sorted(data[label].items()):  # Ahora accedemos a data[label]
                    for dim, f_ref in sorted(dim_values.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
                        dimx, dimy, dimz = dim
                        delta_value = delta.replace('_', '.')  
                        writer.writerow([label,R[label], delta_value, dimx, dimy, dimz, f_ref])

def process_principal_part(output_file, label_struc,R_np, delta_dim_part, aL, k_aL, F):
    structure = os.getcwd()
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("cell, radius [nm],delta,dimx,dimy,dimz,F_value\n")
    delta_list = sorted({entry["delta"] for entry in delta_dim_part if entry["delta"] is not None})
    for delta in delta_list:
        round_value = int(np.round(float(aL/k_aL) / float(delta)))
        dims = []
        dims_sum_bin = [entry["dim"] for entry in delta_dim_part if entry["delta"] == delta][0]
        #dims_sum_bin = [-8,-7,-6,-5,-4,-3,-2,1,0,1,2]
        for sum_dim in dims_sum_bin:
            dims.append(round_value + int(sum_dim))
        delta_folder = str(delta).replace('.','_')

        for j in dims:
            folder_name = f"delta_{delta_folder}_dim_{j}"
            os.chdir(folder_name)
            
            file_path = os.path.join(os.getcwd(), F+".dat")
            
            if os.path.isfile(file_path):
                try:
                    with open(file_path, "r") as file:
                        first_line = file.readline().strip().split()
                        if len(first_line) > 1:
                            value = first_line[1]  # Segunda columna
                            with open(output_file, "a") as out_file:
                                out_file.write(f"{label_struc},{R_np},{delta},{j},{j},{j},{value}\n")
                        else:
                            print(f"Advertencia: No se encontró un valor en {file_path}")
                except Exception as e:
                    print(f"Error al procesar {file_path}: {e}")
            else:
                print(f"Advertencia: Archivo no encontrado en {file_path}")
            os.chdir("..")

def process_reference_part(output_file, base_folder, cell_part, label_struc, F):
    references_file = os.path.join(base_folder, "tot_references.csv")
    contadores = defaultdict(lambda: defaultdict(int))
    delta_map = defaultdict(lambda: defaultdict(set))

    # Leer el archivo de referencias y construir los mapas
    with open(references_file, "r", encoding="utf-8") as ref_file:
        reader = csv.reader(ref_file)
        next(reader)  # Saltar encabezado
        for parts in reader:
            if len(parts) >= 7:
                try:
                    R = parts[0].strip()
                    contador = int(parts[1].strip())
                    pos_tuple = (parts[2], parts[3], parts[4])
                    delta = parts[5].strip().replace(",", ".")
                    dimx = parts[6].strip().replace(",", ".")
                    dimy = parts[7].strip().replace(",", ".")
                    dimz = parts[8].strip().replace(",", ".")
                    key = parts[-1].strip()

                    dim = (dimx,dimy,dimz)
                    # Ahora guardamos el contador asociado a (key, dim)
                    contadores[delta][(key, dim)] = contador
                    delta_map[(label_struc, delta)][key].add(dim)

                except ValueError:
                    print(f"Error de formato en línea -> {','.join(parts)}")

    # Crear archivo de salida si no existe
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("cell,radius [nm],delta,dimx,dimy,dimz,F_reference\n")

    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

    for (label_struc, delta), key_map in delta_map.items():
        delta_name = delta.replace('.', '_')
        delta_folder = os.path.join(base_folder, f"delta_{delta_name}")
        for key, dims in key_map.items():
            key_folder = os.path.join(delta_folder, key)

            # Verificar si el archivo `F_tot_gcanon.dat` existe
            file_path = os.path.join(key_folder, F+".dat")
            if os.path.isfile(file_path):
                try:
                    with open(file_path, "r") as file:
                        value = float(file.readline().strip().split()[1])  # Segunda columna
                        for dim in set(dims):
                            multiplicador = contadores.get(delta, {}).get((key, dim), 1)
                            data[label_struc][delta][dim] += value * multiplicador

                except Exception as e:
                    print(f"Error al procesar {file_path}: {e}")
            else:
                print(f"Advertencia: Archivo no encontrado en {file_path}")

    # Escribir los resultados en el archivo de salida
    for label_struc in cell_part:
        if label_struc in data:  # Verificar si hay datos para esa estructura (fcc o bcc)
            with open(output_file, "a", newline="") as csvfile:
                writer = csv.writer(csvfile)
                for delta, dim_values in sorted(data[label_struc].items()):  # Acceder a los datos de la estructura
                    for dim, f_ref in sorted(dim_values.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
                        dimx, dimy, dimz = dim
                        delta_value = delta.replace('_', '.')  # Asegurarse de que el delta se escribe correctamente
                        writer.writerow([label_struc, R, delta_value, dimx, dimy, dimz, f_ref])  # Escribir en el archivo CSV

