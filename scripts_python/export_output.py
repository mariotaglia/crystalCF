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

def process_principal(output_file, delta_list, aL, F, dims_sum_bin):
    structure = os.getcwd()
    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("delta,dim,F_value\n")

    for delta in delta_list:
        round_value = int(np.round(float(aL) / float(delta)))
        dims = []
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

def process_reference_bin(output_file, dir_fuente, F):
    references_file = os.path.join(dir_fuente, "binary_ref", "tot_references.csv")
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    contadores = defaultdict(lambda: defaultdict(int))
    delta_map = defaultdict(lambda: defaultdict(set))

    # Leer y acumular contadores de las reference
    with open(references_file, "r", encoding="utf-8") as ref_file:
        reader = csv.reader(ref_file)
        next(reader)  # Saltar encabezado
        for parts in reader:
            if len(parts) >= 7:
                try:
                    label = parts[0].strip()  # Primera columna
                    contador = int(parts[1].strip())  
                    pos_tuple = (parts[2], parts[3], parts[4])  # Tupla única de posiciones
                    delta = parts[5].strip().replace(",", ".")  # Quinta columna
                    dim = parts[6].strip().replace(",", ".")  # Sexta columna

                    # Acumular el contador correctamente
                    contadores[delta][(pos_tuple, dim)] = contador  
                    delta_map[(label, delta)][pos_tuple].add(dim)

                except ValueError:
                    print(f"Error de formato en línea -> {','.join(parts)}")

    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("#part,delta,dim,F_reference\n")

    # Procesar los archivos F_tot_gcanon.dat en la estructura existente
    for (label, delta), pos_map in delta_map.items():
        base_folder = os.path.join(dir_fuente, "binary_ref", label)
        delta_name = delta.replace('.', '_')
        delta_folder = os.path.join(base_folder, f"delta_{delta_name}")  

        sub_folder_counter = defaultdict(int)  # Contador para subcarpetas

        for pos_tuple, dims in pos_map.items():
            sorted_dims = sorted(dims, key=lambda x: float(x))  # Asegurar orden numérico
            dim_folder_name = "_".join(f"dim_{dim}" for dim in sorted_dims)
            dim_folder = os.path.join(delta_folder, dim_folder_name)

            sub_folder_counter[dim_folder] += 1  # Contador por carpeta de dimensiones
            sub_folder_name = f"sub_{sub_folder_counter[dim_folder]}"
            sub_folder = os.path.join(dim_folder, sub_folder_name)

            # Verificar si el archivo F_tot_gcanon.dat existe
            file_path = os.path.join(sub_folder, F+".dat")
            if os.path.isfile(file_path):
                try:
                    # Leer y procesar el value de F_tot_gcanon.dat
                    with open(file_path, "r") as file:
                        value = float(file.readline().strip().split()[1])  # Segunda columna

                        # Buscar multiplicador para esta combinación delta/dim
                        for dim in sorted_dims:
                            multiplicador = contadores.get(delta, {}).get((pos_tuple, dim), 1)
                            # Sumar el value multiplicado al diccionario data
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
    # Obtener el nombre de la carpeta actual (que corresponde a la estructura #_ref)
    references_file = os.path.join(base_folder, "tot_references.csv")
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    contadores = defaultdict(lambda: defaultdict(int))
    delta_map = defaultdict(lambda: defaultdict(set))

    # Leer y acumular contadores de las reference
    with open(references_file, "r", encoding="utf-8") as ref_file:
        reader = csv.reader(ref_file)
        next(reader)  # Saltar encabezado
        for parts in reader:
            if len(parts) >= 6:
                try:
                    # Asegurarse de que los datos estén en el formato adecuado
                    contador = int(parts[0].strip())  # Primer columna (entero)
                    pos_tuple = (parts[1], parts[2], parts[3])  # Tupla de posiciones
                    delta = parts[4].strip().replace(",", ".")  # Quinta columna (delta, convertir a formato adecuado)
                    dim = parts[5].strip().replace(",", ".")  # Sexta columna (dim, convertir a formato adecuado)

                    # Acumular el contador correctamente en contadores
                    contadores[delta][(pos_tuple, dim)] = contador  
                    delta_map[(label_struc, delta)][pos_tuple].add(dim)  # Almacenar información de deltas y posiciones

                except ValueError:
                    print(f"Error de formato en línea -> {','.join(parts)}")

    if not os.path.isfile(output_file):
        with open(output_file, "w") as out_file:
            out_file.write("cell,delta,dim,F_reference\n")

    # Procesar los archivos F_tot_gcanon.dat en la estructura existente
    for (label_struc, delta), pos_map in delta_map.items():
        delta_name = delta.replace('.', '_')  # Convertir el delta para usarlo como nombre de carpeta
        delta_folder = os.path.join(base_folder, f"delta_{delta_name}")  

        sub_folder_counter = defaultdict(int)  # Contador para subcarpetas (si existen múltiples subcarpetas)

        # Recorrer cada combinación de posición y dimensiones
        for pos_tuple, dims in pos_map.items():
            sorted_dims = sorted(dims, key=lambda x: float(x))  # Asegurar el orden correcto de las dimensiones
            dim_folder_name = "_".join(f"dim_{dim}" for dim in sorted_dims)  # Nombre de carpeta por dimensiones
            dim_folder = os.path.join(delta_folder, dim_folder_name)

            sub_folder_counter[dim_folder] += 1  # Contar subcarpetas por combinación de dimensiones
            sub_folder_name = f"sub_{sub_folder_counter[dim_folder]}"  # Nombre de la subcarpeta
            
            sub_folder = os.path.join(dim_folder, sub_folder_name)

            # Verificar si el archivo F_tot_gcanon.dat existe en la subcarpeta
            file_path = os.path.join(sub_folder, F+".dat")
            if os.path.isfile(file_path):
                try:
                    # Leer y procesar el valor de F_x.dat
                    with open(file_path, "r") as file:
                        value = float(file.readline().strip().split()[1])  # Segunda columna de F_tot_gcanon.dat

                        # Buscar multiplicador para esta combinación delta/dim
                        for dim in sorted_dims:
                            multiplicador = contadores.get(delta, {}).get((pos_tuple, dim), 1)  # valor por defecto 1
                            # Sumar el valor multiplicado al diccionario data
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

