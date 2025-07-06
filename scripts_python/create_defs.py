import subprocess
import csv
import os
import numpy as np
import shutil
import hashlib
from collections import defaultdict
from transform_refs import calculate_center_ref, process_positions
from function import run_command, read_DEF, write_DEF, extract_references, extract_definitions, generate_references_csv
from itertools import product, combinations
from collections import namedtuple

def transform_reflex(lines, name_bin, structure):
    if name_bin == "Li3Bi":
        n1 = 4
        n2 = 5
    elif name_bin == "NaZn13":
        n1 = 1
        n2 = 32
    elif name_bin == "bccAB6":
        n1 = 2
        n2 = 6

    Particle = namedtuple("Particle", ["center", "semiaxis", "rotation", "coverage", "attraction", "chain_length"])

    # Identificar las secciones
    sections_info = [
        ("! number of particles", 1, 1),
        ("!Center", 1, n1+n2),
        ("!particle semiaxis x y z in nm", 1, n1+n2),
        ("! Rotation", 3, n1+n2),
        ("! coverage", 1, n1+n2),
        ("!Surface-polymer atraction", 1, n1+n2),
        ("!chains lenght", 1, n1+n2)
    ]
    sections_found = {}
    for key, lines_per_particle, tot_particles in sections_info:
        for i, line in enumerate(lines):
            if line.strip() == key:
                sections_found[key] = (i + 1, lines_per_particle, tot_particles)
                break
        else:
            raise ValueError(f"Sección {key} no encontrada")
    new_lines = []
    flag = False; i = 0
    while flag == False: 
        if lines[i].strip() == "! number of particles":
            flag = True
            break
        if lines[i-1].strip() == "!PBC PBC xmin xmax ymin ymax zmin zmax, 1=yes, 2=wall, 0=bulk":
            lines[i] = "PBC 1 1 1 1 1 1\n"
        new_lines.append(lines[i])
        i = i+1

    # Leer partículas
    total = n1 + n2
    particles = []
    for i in range(total):
        center = list(map(float, lines[sections_found["!Center"][0] + i].split()))
        semiaxis = list(map(float, lines[sections_found["!particle semiaxis x y z in nm"][0] + i].split()))
        rotation = [list(map(float, lines[sections_found["! Rotation"][0] + i*3 + j].split())) for j in range(3)]
        coverage = float(lines[sections_found["! coverage"][0] + i])
        attraction = float(lines[sections_found["!Surface-polymer atraction"][0] + i])
        chain_length = int(lines[sections_found["!chains lenght"][0] + i])
        particles.append(Particle(center, semiaxis, rotation, coverage, attraction, chain_length))

    # Transformación
    mirror_planes = [
        (True, False, False),
        (False, True, False),
        (False, False, True),
        (True, True, False),
        (True, False, True),
        (False, True, True),
        (True, True, True),
    ]

    all_particles = []

    for p in particles:
        center_half = [c / 2 for c in p.center]
        all_particles.append(Particle(center_half, p.semiaxis, p.rotation, p.coverage, p.attraction, p.chain_length))

    for plane in mirror_planes:
        for p in particles:
            cx, cy, cz = [c / 2 for c in p.center]
            if plane[0]:
                cx = 1.0 - cx
            if plane[1]:
                cy = 1.0 - cy
            if plane[2]:
                cz = 1.0 - cz
            new_center = [cx, cy, cz]
            all_particles.append(Particle(new_center, p.semiaxis, p.rotation, p.coverage, p.attraction, p.chain_length))

    # Eliminar duplicados
    seen = set()
    unique_particles = []
    for p in all_particles:
        key = tuple(round(x % 1.0, 8) for x in p.center)
        if key not in seen:
            seen.add(key)
            unique_particles.append(p)
    
    #for i, p in enumerate(unique_particles):
    #    print(f"#{i:03d} Pos: {p.center}  Semiejes: {p.semiaxis}  Coverage: {p.coverage}")

    # Guardar archivo
    output_file = os.path.join(structure, "DEFINITIONS.txt")
    with open(output_file, "w") as f:
        f.writelines(new_lines)
        f.write("! number of particles\n")
        f.write(f"{len(unique_particles)}\n")

        f.write("!Center\n")
        for p in unique_particles:
            f.write(f"{p.center[0]:.8f} {p.center[1]:.8f} {p.center[2]:.8f}\n")

        f.write("!particle semiaxis x y z in nm\n")
        for p in unique_particles:
            f.write(f"{p.semiaxis[0]:.8f} {p.semiaxis[1]:.8f} {p.semiaxis[2]:.8f}\n")

        f.write("! Rotation\n")
        for p in unique_particles:
            rot_flat = [f"{val:.6f}" for row in p.rotation for val in row]
            for i in range(0, 9, 3):
                f.write(f"{rot_flat[i]} {rot_flat[i+1]} {rot_flat[i+2]}\n")

        f.write("! coverage\n")
        for p in unique_particles:
            f.write(f"{p.coverage:.2f}\n")

        f.write("!Surface-polymer atraction\n")
        for p in unique_particles:
            f.write(f"{p.attraction:.2f}\n")

        f.write("!chains lenght\n")
        for p in unique_particles:
            f.write(f"{p.chain_length}\n")

def process_principal_binario(reference_DEF, name_bin, delta_dim_bin, aL, n_k_bin, tosubmit, dir_fuente, k_aL, gamma, script_folder):
    structure = os.getcwd()
    pairwise_folder = os.path.join(script_folder,"pairwise_additive_F")
    pairwise_file = "fitpairL12.dat"

    DEF2 =  os.path.join(structure, "DEFINITIONS.txt")
    lines = read_DEF(DEF2)
    #transform_reflex(lines, name_bin, structure)
    DEF =  os.path.join(structure, "DEFINITIONS.txt")
    lines = read_DEF(DEF)
    delta_list = sorted({entry["delta"] for entry in delta_dim_bin if entry["delta"] is not None})
    k = 1
    if name_bin == 'MgZn2':
        k = 2

    for delta in delta_list:
        round_value = int(np.round(float(aL/k_aL["kx"]) / float(delta)))
        dims = []
        dims_sum_bin = [entry["dim"] for entry in delta_dim_bin if entry["delta"] == delta][0]
        #dims_sum_bin = [-8,-7,-6,-5,-4,-3,-2,1,0,1,2]
        for sum_dim in dims_sum_bin:
            dims.append((round_value + int(sum_dim)))
        delta_folder = str(delta).replace('.','_')
        for dim in dims:
            folder_name = f'delta_{delta_folder}_dim_{dim}'
            os.makedirs(folder_name, exist_ok=True)
            os.chdir(folder_name)
            
            shutil.copy(tosubmit, "tosubmit.sh")
            #shutil.copy(f"{pairwise_folder}/{pairwise_file}", pairwise_file)
            shutil.copy(DEF, "DEFINITIONS.txt")
            
            with open("tosubmit.sh", "r") as file:
                content = file.read()
            content = content.replace("_NOMBRE", folder_name)
            with open("tosubmit.sh", "w") as file:
                file.write(content)
            
            with open("DEFINITIONS.txt", "r") as file:
                content = file.read()
            
            #lista = ['Li3Bi', 'NaZn13', 'bccAB6']
            #if name_bin in lista:
            #    content = content.replace("dimx _DIM_", f'dimx {str(dim*2)}').replace("dimy _DIM_", f'dimy {str(int(dim*2*k_aL["kx"]/k_aL["ky"]))}')
            #    content = content.replace("dimz _DIM_", f'dimz {str(int(dim*k*2*k_aL["kx"]/k_aL["kz"]))}').replace("_delta_", str(delta))
            #else:
            content = content.replace("dimx _DIM_", f'dimx {str(dim)}').replace("dimy _DIM_", f'dimy {str(int(dim*k_aL["kx"]/k_aL["ky"]))}')
            content = content.replace("dimz _DIM_", f'dimz {str(int(dim*k*k_aL["kx"]/k_aL["kz"]))}').replace("_delta_", str(delta))

            write_DEF("DEFINITIONS.txt", content)
            lines = read_DEF("DEFINITIONS.txt")
            dir_origen = os.path.abspath(os.path.join(dir_fuente, os.pardir))
            output_DEF_ref = {"part1": os.path.join(dir_origen,"sim_part1","binary_ref","part1"),"part2": os.path.join(dir_fuente,"binary_ref","part2")}
            process_secundario_binario(lines, name_bin, output_DEF_ref, int(dim*k_aL["kx"]), n_k_bin, dir_fuente, delta_list, k_aL, gamma)
            os.chdir(dir_fuente)
            os.chdir(structure)

def process_secundario_binario(lines, name_bin, output_folder, dim, n_k_bin, dir_fuente, delta_bin, k_aL, gamma):
    n1_k_bin = n_k_bin["part1"]; n2_k_bin = n_k_bin["part2"]
    sections_info = [
        ("! number of particles", 1, 1),
        ("!Center", 1, n1_k_bin+n2_k_bin),
        ("!particle semiaxis x y z in nm", 1, n1_k_bin+n2_k_bin),
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
    chain_lenght = {"part1": None, "part2": None}
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
                if key == "!chains lenght":
                    chain_lenght[label] = [elem.replace('\n', '') for elem in new_block][0]
                # Eliminates original params an replace with the new ones
                del modified_lines[start_index: start_index + block_length]
                for line in reversed(new_block):
                    modified_lines.insert(start_index, line)

        output_DEF = os.path.join(output_folder[label], f"DEFINITIONS.txt")
        with open(output_DEF, "w") as f:
            f.writelines(modified_lines)

    for label in ["part1", "part2"]:
        os.chdir(output_folder[label])
        # Extraer datos y procesar posiciones
        data = extract_definitions("DEFINITIONS.txt")
        dimx = int(data.get("dimx"))
        dimy = int(data.get("dimy"))
        dimz = int(data.get("dimz"))
        delta = float(data.get("delta"))
        cdiva = float(data.get("cdiva"))
        centers = data.get("centers", [])
        PBC = data.get("PBC", [])
        R = float(data.get("R")[0])
        nseg = int(chain_lenght[label])
        lseg = float(data["lseg"])
        delta_min = np.min(delta_bin)
        dist_min = (2*R+2*lseg*nseg)
        N_ref = np.round(dist_min*1.50/delta_min)
        N_ref = int(N_ref)
        
        def N_round(N_ref):
            if N_ref%2 == 0:
                N_ref += 1
                return N_ref
            else:
                return N_ref
                
        if k_aL["kx"] == 2:
            Nx = N_round(int(N_ref/1.5))
        else:
            Nx = N_round(int(N_ref*1.2))

        if k_aL["ky"] == 2:
            Ny = N_round(int(N_ref/1.5))
        else:
            Ny = N_round(int(N_ref*1.2))

        if k_aL["kz"] == 2:
            Nz = N_round(int(N_ref/1.5))
        else:
            Nz = N_round(int(N_ref*1.2))

        center_ref_list = calculate_center_ref(N_ref, centers, dimx, dimy, dimz, delta, cdiva, PBC)
        pos_out, _ = process_positions(center_ref_list)

        references = extract_references("references.csv")
        os.chdir(os.path.abspath(os.path.join(output_folder[label], os.pardir)))

        def_ref_path = os.path.join(label, "DEFINITIONS.txt")

        with open(def_ref_path, "r") as file:
            lines = file.readlines()
    
        new_lines = []
        for line in lines:
            if line.startswith("dimx"):
                new_lines.append(f"dimx {int(Nx)}\n")
            elif line.startswith("dimy"):
                new_lines.append(f"dimy {int(Ny)}\n")

            elif line.startswith("dimz"):
                k = 1
                if name_bin == 'MgZn2':
                    k = 2

                new_lines.append(f"dimz {int(Nz*k)}\n")
            elif line.startswith("delta"):
                new_lines.append("delta _delta_\n")
            else:
                new_lines.append(line)
        
        with open(def_ref_path, "w") as file:
            file.writelines(new_lines)
        generate_references_csv(references, os.getcwd(), delta, R, dimx, dimy, dimz, label, name_bin, gamma)

def process_terciario_binario(output_folder, name_bin, references, tosubmit, dir_fuente, n_k_bin):
    '''
    From the references file, separates each particle and create a folder delta/dim/sub
    for each equivalent relative position to estimate the energy of a particle isolated.
    '''
    delta_map = defaultdict(lambda: defaultdict(set))
    for reference in references:
        label = reference[0]
        pos_tuple = (reference[3], reference[4], reference[5])  # Tupla única de posiciones
        delta = str(reference[6]).replace('.', '_')  # Formato correcto de delta
        dim = reference[7]  # Dimensión
        key = reference[10] 
        delta_map[(label, delta)][pos_tuple].add(key)

    # Segunda fase: creación de carpetas y archivos
    for (label, delta), pos_map in delta_map.items():
        DEF_ref = os.path.join(dir_fuente, "binary_ref", label, "DEFINITIONS.txt")
        delta_folder = os.path.join(output_folder, label, f"delta_{delta}")
        os.makedirs(delta_folder, exist_ok=True)

        for pos_tuple,key_value in pos_map.items():
            sorted_key_value = sorted(key_value, key=lambda x: str(x))  
            key_folder_name = "_".join(f"{key_value}" for key_value in sorted_key_value)
            key_folder = os.path.join(delta_folder, key_folder_name)
            os.makedirs(key_folder, exist_ok=True)
            
            # Copiar y modificar tosubmit.sh
            shutil.copy(tosubmit, os.path.join(key_folder, "tosubmit.sh"))
            with open(os.path.join(key_folder, "tosubmit.sh"), "r") as file:
                content = file.read()
            content = content.replace("_NOMBRE", key_folder)
            with open(os.path.join(key_folder, "tosubmit.sh"), "w") as file:
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
            k = 1
            if name_bin == 'MgZn2':
                k = 2
            content = content.replace("dimx _DIM_", f'dimx {str(N_ref)}').replace("dimy _DIM_", f'dimy {str(N_ref)}')
            content = content.replace("dimz _DIM_", f'dimz {str(N_ref*k)}').replace("_delta_", str(delta_num))

            with open(os.path.join(key_folder, "DEFINITIONS.txt"), "w") as file:
                file.write(content)
            
            lines = read_DEF(os.path.join(key_folder, "DEFINITIONS.txt"))
            definitions_ref_edit(lines, key_folder, *pos_tuple, n_k_bin[label])

def definitions_ref_edit(lines, ref_folder, pos1, pos2, pos3, ni_k_bin):
    pos1, pos2, pos3 = pos1, pos2, pos3
    sections_info = {
        "! number of particles": 1,
        "!save vtk?": 1,
        "!Center": ni_k_bin,  
        "!particle semiaxis x y z in nm": ni_k_bin,
        "! Rotation": 3*ni_k_bin,
        "! coverage": ni_k_bin,
        "!Surface-polymer atraction": ni_k_bin,
        "!chains lenght": ni_k_bin
    }
    
    modified_lines = lines.copy()
    for i, line in enumerate(modified_lines):
        if line.strip() == "! number of particles":
            modified_lines[i+1] = "1\n"
            break
        if line.strip() == "!save vtk?":
            modified_lines[i+1] = "vtkflag 0\n"
            continue
    
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

