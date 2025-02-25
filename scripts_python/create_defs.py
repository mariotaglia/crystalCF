import subprocess
import csv
import os
import numpy as np
import shutil
from collections import defaultdict
from transform_refs import calculate_center_ref, process_positions
from function import run_command, read_DEF, write_DEF, extract_definitions, extract_references

def process_principal_binario(reference_DEF, delta_list, aL, n, k_bin, tosubmit, dir_fuente, dims_sum_bin):
    structure = os.getcwd()
    DEF =  os.path.join(structure, "DEFINITIONS.txt")
    lines = read_DEF(DEF)

    for delta in delta_list:
        round_value = int(np.round(float(aL) / float(delta)))
        dims = []
        for sum_dim in dims_sum_bin:
            dims.append(round_value + int(sum_dim))
        delta_folder = str(delta).replace('.','_')
        for dim in dims:
            folder_name = f"delta_{delta_folder}_dim_{dim}"
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
            content = content.replace("_DIM_", str(dim)).replace("_delta_", str(delta))
            write_DEF("DEFINITIONS.txt", content)
            lines = read_DEF("DEFINITIONS.txt")

            output_DEF_ref = os.path.join(dir_fuente,"binary_ref")
            process_secundario_binario(lines, output_DEF_ref, delta, dim, n, k_bin, dir_fuente, delta_list)
            os.chdir(dir_fuente)
            os.chdir(structure)

def process_secundario_binario(lines, output_folder, delta, dim, n, k_bin, dir_fuente, delta_bin):
    n1 = n["part1"]; n2 = n["part2"]
    sections_info = [
        ("! number of particles", 1, 1),
        ("!Center", 1, k_bin*(n1+n2)),
        ("!particle semiaxis x y z in nm", 1, k_bin*(n1+n2)),
        ("! Rotation", 3, k_bin*(n1+n2)),
        ("! coverage", 1, k_bin*(n1+n2)),
        ("!Surface-polymer atraction", 1, k_bin*(n1+n2))
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

    configs = [("part1", 0, k_bin*n1), ("part2", k_bin*n1, k_bin*n2)]  # (Name, Offset, number of particles)

    for label, offset, num_particles in configs:
        modified_lines = lines.copy()

        for key, start_index, lines_per_particle, tot_particles in sorted(sections_found, key=lambda x: x[1], reverse=True):
            if key == "! number of particles":
                modified_lines[start_index] = f"{num_particles}\n"
            else:
                block_length = tot_particles * lines_per_particle
                start_offset = start_index + offset*lines_per_particle
                end_offset = start_offset + num_particles*lines_per_particle
                # Extract the particle params from DEFINITIONS
                new_block = modified_lines[start_offset:end_offset]

                # Eliminates original params an replace with the new ones
                del modified_lines[start_index: start_index + block_length]
                for line in reversed(new_block):
                    modified_lines.insert(start_index, line)

        output_DEF = os.path.join(output_folder, f"{label}/DEFINITIONS.txt")
        with open(output_DEF, "w") as f:
            f.writelines(modified_lines)

    for label in ["part1", "part2"]:
        os.chdir(os.path.join(dir_fuente, "binary_ref",label))

        # Extraer datos y procesar posiciones
        data = extract_definitions("DEFINITIONS.txt")
        dimx = int(data.get("dimx"))
        dimy = int(data.get("dimy"))
        dimz = int(data.get("dimz"))
        delta = float(data.get("delta"))
        cdiva = float(data.get("cdiva"))
        centers = data.get("centers", [])
        R = float(data.get("R")[0])
        nseg = int(data.get("nseg")); lseg = float(data.get("lseg"))
        delta_min = np.min(delta_bin)
        N_ref = np.round((2*R+2*lseg*nseg)*1.50/delta_min)
        N_ref = int(N_ref)
        if N_ref%2 == 0:
            N_ref += 1

        center_ref_list = calculate_center_ref(N_ref, centers, dimx, dimy, dimz, delta, cdiva)
        pos_out, _ = process_positions(center_ref_list)

        references = extract_references("references.csv")
        os.chdir(os.path.join(dir_fuente, "binary_ref"))

        def_ref_path = os.path.join(dir_fuente, "binary_ref",label, "DEFINITIONS.txt")

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
        generate_references_csv(references, os.getcwd(), delta, dim, label)

def process_terciario_binario(output_folder, references, tosubmit, dir_fuente, n, k_bin):
    '''
    From the references file, separates each particle and create a folder delta/dim/sub
    for each equivalent relative position to estimate the energy of a particle isolated.
    '''
    delta_map = defaultdict(lambda: defaultdict(set)) 
    for reference in references:
        label = reference[0]
        pos_tuple = (reference[2], reference[3], reference[4])  # Tupla única de posiciones
        delta = str(reference[5]).replace('.', '_')  # Formato correcto de delta
        dim = reference[6]  # Dimensión
        delta_map[(label, delta)][pos_tuple].add(dim)

    # Segunda fase: creación de carpetas y archivos
    for (label, delta), pos_map in delta_map.items():
        DEF_ref = os.path.join(dir_fuente, "binary_ref", label, "DEFINITIONS.txt")
        delta_folder = os.path.join(output_folder, label, f"delta_{delta}")  # Nombre correcto
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
            content = content.replace("_NOMBRE", sub_folder_name)
            with open(os.path.join(sub_folder, "tosubmit.sh"), "w") as file:
                file.write(content)

            if not os.path.exists(DEF_ref):
                print(f"ERROR: {DEF_ref} not found.")
                continue
            
            delta_num = delta.replace('_', '.')
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
            
            lines = read_DEF(os.path.join(sub_folder, "DEFINITIONS.txt"))
            definitions_ref_edit(lines, sub_folder, *pos_tuple, n[label]*k_bin)

def generate_references_csv(references, output_folder, delta_value, dim_value, label):
    references_path = os.path.join(output_folder, "tot_references.csv")
    if not os.path.exists(references_path):
        with open(references_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['#part'] + references[0] + ['delta', 'dim'])
            for row in references[1:]:
                writer.writerow([label] + row + [delta_value, dim_value])
    else:
        with open(references_path, mode='a', newline='') as file:
            writer = csv.writer(file)
            for row in references[1:]:
                writer.writerow([label] + row + [delta_value, dim_value])
    return references_path

def definitions_ref_edit(lines, ref_folder, pos1, pos2, pos3, ni_k_bin):
    pos1, pos2, pos3 = pos1, pos2, pos3
    sections_info = {
        "! number of particles": 1,
        "!Center": ni_k_bin,  
        "!particle semiaxis x y z in nm": ni_k_bin,
        "! Rotation": 3*ni_k_bin,
        "! coverage": ni_k_bin,
        "!Surface-polymer atraction": ni_k_bin
    }
    
    modified_lines = lines.copy()
    for i, line in enumerate(modified_lines):
        if line.strip() == "! number of particles":
            modified_lines[i+1] = "1\n"
            break
    
    for key, num_lines in sections_info.items():
        for i, line in enumerate(modified_lines):
            if line.strip() == key:
                modified_lines[i+1:i+1 + ((num_lines // ni_k_bin) * (ni_k_bin - 1))] = []
                break
    
    for i, line in enumerate(modified_lines):
        if line.strip() == "!Center":
            modified_lines[i+1] = f"{pos1} {pos2} {pos3}\n"
            break
    
    output_DEF = os.path.join(ref_folder, "DEFINITIONS.txt")
    with open(output_DEF, "w") as f:
        f.writelines(modified_lines) 
