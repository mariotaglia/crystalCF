import os
import shutil
import csv
import numpy as np
import subprocess
from collections import defaultdict
from transform_refs import extract_definitions, calculate_center_ref, process_positions

def run_command(command):
    """Run the command on the terminal and return it's output."""
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    return result.stdout.strip()

def read_DEF(file_path):
    """Extract the lines from DEF."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines

def write_DEF(file_path, lines):
    """Export DEF with lines."""
    with open(file_path, 'w') as f:
        f.writelines(lines)

def extract_R_part(definitions_path):
    """Extract the first radius from DEF"""
    R_np = None
    lines = read_DEF(definitions_path)
    size_index = None
    for i, line in enumerate(lines):
        if line.strip() == "!particle semiaxis x y z in nm":
            size_index = i + 1
            break

    if size_index is not None:
        try:
            R_np = float(lines[size_index].split()[0])
        except (ValueError, IndexError):
            print("Error al leer R de la part.")
    return R_np

def extract_references(csv_path):
    """Extract data from csv."""
    references = []
    with open(csv_path, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            references.append(row)
    return references

def generate_references_part_csv(references, output_folder, delta_value, dim_value):
    references_path = os.path.join(output_folder, "tot_references.csv")
    # Verifica si el archivo ya existe, si no lo crea
    if not os.path.exists(references_path):
        with open(references_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(references[0] + ['delta', 'dim'])  # Agregar encabezado
            for row in references[1:]:
                writer.writerow(row + [delta_value, dim_value])
    else:
        with open(references_path, mode='a', newline='') as file:
            writer = csv.writer(file)
            for row in references[1:]:
                writer.writerow(row + [delta_value, dim_value])

def path_carpeta(folder_init,n):
    x = os.path.dirname(folder_init)
    for i in range(n-1):
        x = os.path.dirname(x)
    return x

def process_principal_part(reference_DEF, delta_list, aL, tosubmit, dir_fuente):
    structure = os.getcwd()
    DEF =  os.path.join(structure, "DEFINITIONS.txt")
    lines = read_DEF(DEF)

    for delta in delta_list:
        round_value = int(np.round(float(aL) / float(delta)))
        dims = [round_value - 1, round_value, round_value + 1]
        delta_folder = str(delta).replace('.','_')
        for j in dims:
            folder_name = f"delta_{delta_folder}_dim_{j}"
            os.makedirs(folder_name, exist_ok=True)
            os.chdir(folder_name)
            
            shutil.copy(tosubmit, "tosubmit.sh")
            shutil.copy(DEF, "DEFINITIONS.txt")
            
            with open("tosubmit.sh", "r") as file:
                content = file.read()
            content = content.replace("_NOMBRE", folder_name)
            with open("tosubmit.sh", "w") as file:
                file.write(content)
            
            with open("DEFINITIONS.txt", "r") as file:
                content = file.read()
            content = content.replace("_DIM_", str(j)).replace("_delta_", str(delta))
            write_DEF("DEFINITIONS.txt", content)
            lines = read_DEF("DEFINITIONS.txt")

            output_DEF_ref = os.path.join(dir_fuente,"binary_ref")
            process_secundario_part(structure, j, dir_fuente, delta_list)
            os.chdir(dir_fuente)
            os.chdir(structure)

def process_secundario_part(struc, dim, dir_fuente_part, delta_list):
    data = extract_definitions("DEFINITIONS.txt")
    dimx = int(data.get("dimx"))
    dimy = int(data.get("dimy"))
    dimz = int(data.get("dimz"))
    delta = float(data.get("delta"))
    cdiva = float(data.get("cdiva"))
    centers = data.get("centers", [])
    R = float(data.get("R")[0])
    nseg = int(data.get("nseg")); lseg = float(data.get("lseg"))
    delta_min = np.min(delta_list)

    N_ref = np.round((2*R+2*lseg*nseg)*1.50/delta_min)
    N_ref = int(N_ref)
    if N_ref%2 == 0:
        N_ref += 1

    center_ref_list = calculate_center_ref(N_ref, centers, dimx, dimy, dimz, delta, cdiva)
    pos_out, _ = process_positions(center_ref_list)

    references = extract_references("references.csv")
    def_ref_path = os.path.join(dir_fuente_part, f"{struc}_ref", "DEFINITIONS.txt")
    shutil.copy("DEFINITIONS.txt", def_ref_path)
    with open(def_ref_path, "r") as file:
        lines = file.readlines()

    new_lines = []
    for line in lines:
        if line.startswith("dimx"):
            new_lines.append(f"dimx {N_ref}\n")
        elif line.startswith("dimy"):
            new_lines.append(f"dimy {N_ref}\n")
        elif line.startswith("dimz"):
            new_lines.append(f"dimz {N_ref}\n")
        elif line.startswith("delta"):
            new_lines.append("delta _delta_\n")
        else:
            new_lines.append(line)

    with open(def_ref_path, "w") as file:
        file.writelines(new_lines)
    os.chdir(os.path.join(dir_fuente_part, f"{struc}_ref"))
    generate_references_part_csv(references, os.getcwd(), delta, dim)

def process_terciario_part(output_folder, references, DEF, tosubmit):
    delta_map = defaultdict(lambda: defaultdict(set))
    DEF_ref = os.path.join(output_folder, "DEFINITIONS.txt")

    # Agrupar referencias por delta, dimensiones y posiciones únicas
    for reference in references:
        pos_tuple = (reference[1], reference[2], reference[3])  # Tupla única de posiciones
        delta = str(reference[4]).replace('.', '_')  # Formato correcto de delta
        dim = reference[5]  # Dimensión
        delta_map[delta][pos_tuple].add(dim)

    # Segunda fase: creación de carpetas y archivos
    for delta, pos_map in delta_map.items():
        delta_folder = os.path.join(output_folder, f"delta_{delta}")  
        os.makedirs(delta_folder, exist_ok=True)
        sub_folder_counter = defaultdict(int)  # Contador para subcarpetas

        for pos_tuple, dims in pos_map.items():
            sorted_dims = sorted(dims, key=lambda x: float(x))  # Asegurar orden numérico 
            dim_folder_name = "_".join(f"dim_{dim}" for dim in sorted_dims)
            dim_folder = os.path.join(delta_folder, dim_folder_name)
            os.makedirs(dim_folder, exist_ok=True)
            
            sub_folder_counter[dim_folder] += 1  # Contador por carpeta de dimensiones
            sub_folder_name = f"sub_{sub_folder_counter[dim_folder]}"
            sub_folder = os.path.join(dim_folder, sub_folder_name)
            os.makedirs(sub_folder, exist_ok=True)

            # Copiar y modificar tosubmit.sh

            shutil.copy(tosubmit, os.path.join(sub_folder, "tosubmit.sh"))
            with open(os.path.join(sub_folder, "tosubmit.sh"), "r") as file:
                content = file.read()
            content = content.replace("_NOMBRE", sub_folder)
            with open(os.path.join(sub_folder, "tosubmit.sh"), "w") as file:
                file.write(content)

            # Modificar y escribir DEFINITIONS.txt
            delta_num = delta.replace('_', '.')  # Restaurar el formato de delta
            lines = read_DEF(DEF_ref)
            for i, line in enumerate(lines):
                line = lines[i].strip()
                parts = line.split()
                if "dimx" in line:
                    N_ref = str(int(parts[1]))

            with open(DEF_ref, "r") as file:
                content = file.read()
            content = content.replace("_DIM_", N_ref).replace("_delta_", str(delta_num))
            with open(os.path.join(sub_folder, "DEFINITIONS.txt"), "w") as file:
                file.write(content)
            
            # Editar DEFINITIONS.txt con referencias de posiciones
            with open(os.path.join(sub_folder, "DEFINITIONS.txt"), "r") as file:
                lines = file.readlines()

            pos1, pos2, pos3 = pos_tuple
            definitions_ref_edit(lines, sub_folder, pos1, pos2, pos3)

def generate_references_part_csv(references, output_folder, delta_value, dim_value):
    references_path = os.path.join(output_folder, "tot_references.csv")
    # Verifica si el archivo ya existe, si no lo crea
    if not os.path.exists(references_path):
        with open(references_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(references[0] + ['delta', 'dim'])  # Agregar encabezado
            for row in references[1:]:
                writer.writerow(row + [delta_value, dim_value])
    else:
        with open(references_path, mode='a', newline='') as file:
            writer = csv.writer(file)
            for row in references[1:]:
                writer.writerow(row + [delta_value, dim_value])

def definitions_ref_edit(lines, ref_folder, pos1, pos2, pos3):
    pos1, pos2, pos3 = pos1, pos2, pos3
    sections_info = {
        "! number of particles": 1,  # Solo 1 partícula
        "!Center": 2,  
        "!particle semiaxis x y z in nm": 2,
        "! Rotation": 6,  # 3 líneas por partícula, 2 partículas en original
        "! coverage": 2,
        "!Surface-polymer atraction": 2
    }
    
    # Encontrar las posiciones de cada sección en el archivo
    sections_found = {}
    for key, num_lines in sections_info.items():
        for i, line in enumerate(lines):
            if line.strip() == key:
                sections_found[key] = i + 1  # Guardamos la posición inicial de la sección
                break 

    modified_lines = lines.copy()
    for i, line in enumerate(modified_lines):
        if line.strip() == "! number of particles":
            check_num_part = int(lines[i+1].strip())
            modified_lines[i+1] = "1\n"
            break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "!Center":
            # Reemplazar la primera línea con la posición de la partícula (los valores pueden ser flotantes)
            modified_lines[i+1] = f"{pos1} {pos2} {pos3}\n"
            # Eliminar la segunda línea de centro (si existe y no es un encabezado)
            if i+2 < len(modified_lines) and not modified_lines[i+2].strip().startswith("!"):
                del modified_lines[i+2]
            break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "!particle semiaxis x y z in nm":
            if i+2 < len(modified_lines) and not modified_lines[i+2].strip().startswith("!"):
                del modified_lines[i+2]
            break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "! Rotation":
            # Se asume que las siguientes 6 líneas corresponden a dos rotaciones (3 líneas cada una)
            # Eliminamos las líneas correspondientes a la segunda partícula.
            if i+6 < len(modified_lines):
                # Las líneas a conservar son i+1, i+2 y i+3
                if check_num_part>1:
                    del modified_lines[i+4 : i+7]  # elimina líneas en posiciones i+4, i+5 y i+6
            break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "! coverage":
            if i+2 < len(modified_lines) and not modified_lines[i+2].strip().startswith("!"):
                del modified_lines[i+2]
            break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "!Surface-polymer atraction":
            if i+2 < len(modified_lines) and not modified_lines[i+2].strip().startswith("!"):
                del modified_lines[i+2]
            break
    
    # Guardar el archivo DEFINITIONS.txt en la carpeta correspondiente
    output_DEF = os.path.join(ref_folder, "DEFINITIONS.txt")
    with open(output_DEF, "w") as f:
        f.writelines(modified_lines) 
