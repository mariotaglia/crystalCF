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

def process_principal(output_file, delta_dim_bin, aL, F):
    structure = os.getcwd()
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("delta,dim,F_value\n")

    delta_list = sorted({entry["delta"] for entry in delta_dim_bin if entry["delta"] is not None})
    for delta in delta_list:
        round_value = int(np.round(float(aL) / float(delta)))
        dims = []
        dims_sum_bin = [entry["dim"] for entry in delta_dim_bin if entry["delta"] == delta][0]
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
                            value = first_line[1]  # Segunda columna
                            with open(output_file, "a") as out_file:
                                out_file.write(f"{delta},{dim},{value}\n")
                        else:
                            print(f"Advertencia: No se encontró un valor en {file_path}")
                except Exception as e:
                    print(f"Error al procesar {file_path}: {e}")
            else:
                print(f"Advertencia: Archivo no encontrado en {file_path}")
            os.chdir("..")

def process_reference_bin(output_file, dir_inicial, F):
    dir_fuente = {"part1": os.path.join(dir_inicial,"sim_part1"),"part2": os.path.join(os.getcwd())}
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

    for label in ["part1", "part2"]:
        references_file = os.path.join(dir_fuente[label], "binary_ref", "tot_references.csv")
        contadores = defaultdict(lambda: defaultdict(int))
        delta_map = defaultdict(lambda: defaultdict(set))

        # Leer el archivo de referencias y construir los mapas
        with open(references_file, "r", encoding="utf-8") as ref_file:
            reader = csv.reader(ref_file)
            next(reader)  # Saltar encabezado
            for parts in reader:
                if len(parts) >= 8:
                    try:
                        label = parts[0].strip()
                        contador = int(parts[1].strip())
                        pos_tuple = (parts[2], parts[3], parts[4])
                        delta = parts[5].strip().replace(",", ".")
                        dim = parts[6].strip().replace(",", ".")
                        key = parts[7].strip()

                        # Ahora guardamos el contador asociado a (key, dim)
                        contadores[delta][(key, dim)] = contador
                        delta_map[(label, delta)][key].add(dim)

                    except ValueError:
                        print(f"Error de formato en línea -> {','.join(parts)}")

        # Crear archivo de salida si no existe
        if not os.path.isfile(output_file):
            with open(output_file, "w") as out_file:
                out_file.write("#part,delta,dim,F_reference\n")

        for (label, delta), key_map in delta_map.items():
            base_folder = os.path.join(dir_fuente[label], "binary_ref", label)
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
                            for dim in dims:
                                multiplicador = contadores.get(delta, {}).get((key, dim), 1)
                                data[label][delta][dim] += value * multiplicador

                    except Exception as e:
                        print(f"Error al procesar {file_path}: {e}")
                else:
                    print(f"Advertencia: Archivo no encontrado en {file_path}")

    # Escribir los resultados en el archivo de salida
    for label in ["part1", "part2"]:
        if label in data:  # Verificar si hay datos para esa etiqueta
            with open(output_file, "a", newline="") as csvfile:
                writer = csv.writer(csvfile)
                for delta, dim_values in sorted(data[label].items()):  # Ahora accedemos a data[label]
                    for dim, f_ref in sorted(dim_values.items(), key=lambda x: float(x[0])):
                        delta_value = delta.replace('_', '.')  
                        writer.writerow([label, delta_value, dim, f_ref])

def process_principal_part(output_file, label_struc, delta_list, aL, F):
    structure = os.getcwd()
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("cell,delta,dim,F_value\n")

    for delta in delta_list:
        delta_folder = str(delta).replace('.', '_')
        round_value = int(np.round(float(aL) / float(delta)))
        dims = [round_value - 1, round_value, round_value + 1]
        
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
                                out_file.write(f"{label_struc},{delta},{j},{value}\n")
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
                    contador = int(parts[0].strip())
                    pos_tuple = (parts[1], parts[2], parts[3])
                    delta = parts[4].strip().replace(",", ".")
                    dim = parts[5].strip().replace(",", ".")
                    key = parts[6].strip()

                    # Ahora guardamos el contador asociado a (key, dim)
                    contadores[delta][(key, dim)] = contador
                    delta_map[(label_struc, delta)][key].add(dim)

                except ValueError:
                    print(f"Error de formato en línea -> {','.join(parts)}")

    # Crear archivo de salida si no existe
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("cell,delta,dim,F_reference\n")

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
                        for dim in dims:
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
                    for dim, f_ref in sorted(dim_values.items(), key=lambda x: float(x[0])):
                        delta_value = delta.replace('_', '.')  # Asegurarse de que el delta se escribe correctamente
                        writer.writerow([label_struc, delta_value, dim, f_ref])  # Escribir en el archivo CSV

