import subprocess
import csv
import os
import shutil
import pandas as pd
import glob
import re
import numpy as np
import hashlib
from references.dependecies_init import edit_def_bin, extract_definitions
from collections import defaultdict
import math

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
        
def extract_params_init(params_init, cond):
    data = {
        "name": None, "n1": None, "n2": None, "R1": None,
        "gamma list": [], "list delta bin": [], "list gamma delta sum dim": [],
        "cell part": [], "list delta part": {}, "list part delta sum dim": [], "num cell part": {}, "cell bin factor": None,
        "num cell bin": None, "flag generate energy vs aL curves": None,
        "flag reflexion binary": None, "flag reflexion part": None, "PBC": []
    }

    lines = read_DEF(params_init)
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        if line == "!name":
            data["name"] = lines[i + 1].strip("\n")
            i += 1
        elif line.startswith(("n1", "n2")):
            key, value = line.split()
            data[key] = int(value)
        elif line == "!radius part1":
            data["R1"] = float(lines[i+1].split()[1])
            i += 1
        elif line == "!aL cell bin factor":
            data["cell bin factor"] = eval(lines[i+1].split()[0])
            i += 1
        elif line == "!list gamma":
            data["gamma list"] = [float(x) for x in lines[i+1].strip("[]\n").split(",")]
            i += 1
        elif line == "!list delta bin":
            data["list delta bin"] = [float(x) for x in lines[i+1].strip("[]\n").split(",")]
            i += 1
        elif line == "!point sum dims bin":
            j = i + 1
            ref = [-1, 0, 1]
            delta_ref = data["list delta bin"]
            gamma_delta_map = {(g, d): dims for g, d, dims in data["list gamma delta sum dim"]}  # Mapeo para evitar duplicados
            
            while j < len(lines) and not lines[j].startswith("!"):
                parts = lines[j].split(maxsplit=2)  # ["gamma", value, list] o ["gamma", value, list, "delta", delta_value]
                if len(parts) >= 3 and parts[0] == "gamma":
                    gamma_value = float(parts[1])
                    if "delta" in parts[2]:
                        sum_dim_part, delta_part = parts[2].split("delta")
                        sum_dim_values = [int(x) for x in re.findall(r'-?\d+', parts[2])[:-2]]
                        delta_value = float(delta_part.strip())
                    else:
                        delta_value = None  # No hay delta explícito
                        sum_dim_values = [int(x) for x in re.findall(r'-?\d+', parts[2])]
                    
                    if delta_value is None:
                        for delta in delta_ref:
                            gamma_delta_map[(gamma_value, delta)] = sum_dim_values  # Reemplazo si existe
                    else:
                        gamma_delta_map[(gamma_value, delta_value)] = sum_dim_values  # Reemplazo si existe
                
                j += 1
            i = j - 1

            for gamma in data["gamma list"]:
                for delta in delta_ref:
                    if (gamma, delta) not in gamma_delta_map:
                        gamma_delta_map[(gamma, delta)] = ref.copy()
            
            # Convertir el diccionario de nuevo en la lista con la estructura adecuada
            data["list gamma delta sum dim"] = [{"gamma": g, "delta": d, "dim": dims} for (g, d), dims in gamma_delta_map.items()]
            data["list gamma delta sum dim"].sort(key=lambda x: (x["gamma"], x["delta"]))

        elif line == "!num cell bin":
            data["num cell bin"] = int(lines[i+1].split()[1])
            i += 1
        
        elif line == "!cell part":
            j = i + 1
            while j < len(lines) and not lines[j].startswith("!"):
                cell_name = lines[j].strip()
                if cell_name:  # Verifica que la línea no esté vacía
                    data["cell part"].append(cell_name)
                j += 1
            i = j - 1  # Saltar líneas procesadas
        
        elif line == "!list delta part":
            j = i + 1
            while j < len(lines) and not lines[j].startswith("!"):
                parts = lines[j].split(maxsplit=1)
                if len(parts) > 1 and parts[0] in data["cell part"]:
                    values = parts[1].strip("[]\n").replace(",", " ").split()
                    data["list delta part"][parts[0]] = [float(x) for x in values]
                j += 1
            i = j - 1

        elif line == "!point sum dims part":
            j = i + 1
            ref = [-1, 0, 1]
            struct_ref = ["fcc", "bcc"]
            delta_ref = data["list delta part"]
            part_delta_map = {(p, s, d): dims for p, s, d, dims in data["list part delta sum dim"]}

            while j < len(lines) and not lines[j].startswith("!"):
                line_j = lines[j].strip()
                if line_j.startswith("part"):
                    parts = line_j.split(maxsplit=3)
                    part_value = parts[1]
                    part = "part1" if part_value == "A" else "part2"

                    if len(parts) < 3:
                        j += 1
                        continue  # Línea incompleta

                    token_2 = parts[2]
                    if token_2 in struct_ref:
                        cell_values = [token_2]
                        raw_info = parts[3] if len(parts) > 3 else ""
                    else:
                        cell_values = struct_ref.copy()
                        raw_info = line_j.split(maxsplit=2)[2]  # todo lo que viene después de 'part A'

                    # Extraer delta y lista
                    if "delta" in raw_info:
                        split_part = raw_info.split("delta")
                        dims_text = split_part[0]
                        delta_text = split_part[1]
                        sum_dim_values = [int(x) for x in re.findall(r'-?\d+', dims_text)]
                        try:
                            delta_value = float(delta_text.strip())
                        except ValueError:
                            raise ValueError(f"No se puede interpretar delta en la línea: {line_j}")
                    else:
                        sum_dim_values = [int(x) for x in re.findall(r'-?\d+', raw_info)]
                        delta_value = None

                    # Asignar
                    for cell in cell_values:
                        if delta_value is None:
                            for delta in delta_ref[cell]:
                                part_delta_map[(part, cell, delta)] = sum_dim_values
                        else:
                            part_delta_map[(part, cell, delta_value)] = sum_dim_values
                j += 1
            i = j - 1

            # Completar con valores por defecto
            for part in ["part1", "part2"]:
                for cell in struct_ref:
                    for delta in delta_ref[cell]:
                        if (part, cell, delta) not in part_delta_map:
                            part_delta_map[(part, cell, delta)] = ref.copy()

            # Guardar en formato lista
            data["list part delta sum dim"] = [
                {"part": p, "cell": c, "delta": d, "dim": dims}
                for (p, c, d), dims in part_delta_map.items()
            ]
            data["list part delta sum dim"].sort(key=lambda x: (x["part"], x["cell"], x["delta"]))

        elif line == "!num cell part":
            j = i + 1
            while j < len(lines) and not lines[j].startswith("!"):
                parts = lines[j].split()
                if len(parts) == 2 and parts[0].startswith("k_"):
                    key = parts[0][2:]  # Deletes "k_"
                    if key in data["cell part"]:
                        data["num cell part"][key] = int(parts[1])
                j += 1
            i = j - 1  # Saltar líneas procesadas
        elif line == "!flag generate energy vs aL curves":
            value = lines[i + 1].strip("\n")
            if value == "True":
                data["flag generate energy vs aL curves"] = True
            else:
                data["flag generate energy vs aL curves"] = False
            i += 1
        elif line == "!flag use reflexion planes binary":
            value = lines[i + 1].strip("\n")
            if value == "True":
                data["flag reflexion binary"] = True
            else:
                data["flag reflexion binary"] = False
            i += 1
        elif line == "!flag use reflexion planes part":
            value = lines[i + 1].strip("\n")
            if value == "True":
                data["flag reflexion part"] = True
            else:
                data["flag reflexion part"] = False
            i += 1

        i += 1

    if data["flag reflexion binary"] == True and cond == False:
        if data["name"] == "MgZn2":
            data["n1"] = 2
            data["n2"] = 7
            data["num cell bin"] = 1
        elif data["name"] == "NaCl" or data["name"] == "CsCl":
            data["n1"] = 2
            data["n2"] = 2
            data["num cell bin"] = 1
        elif data["name"] == "Li3Bi":
            data["n1"] = 4
            data["n2"] = 5
            data["num cell bin"] = 1
        elif data["name"] == "NaZn13":
            data["n1"] = 1
            data["n2"] = 32
            data["num cell bin"] = 1
        elif data["name"] == "bccAB6":
            data["n1"] = 2
            data["n2"] = 6
            data["num cell bin"] = 1
    return data

def extract_definitions(definitions_path):
    data = {
        "dimx": None,
        "dimy": None,
        "dimz": None,
        "delta": None,
        "cdiva": None,
        "num_particles": None,
        "centers": [],
        "R": [],
        "lseg": None,
        "PBC": []
    }
    lines = read_DEF(definitions_path)
    for i, line in enumerate(lines):
        line = lines[i].strip()
        parts = line.split()
        if len(parts) == 2 and parts[0] in data:
            try:
                data[parts[0]] = float(parts[1])
            except ValueError:
                print(f"Error al leer {parts[0]}: {parts[1]}")

        elif line == "! number of particles":
            try:
                data["num_particles"] = float(lines[i+1].split()[0])
            except ValueError:
                print("Error al leer el numero de particulas.")
        elif line == "!cdiva":
            try:
                data["cdiva"] = float(lines[i+1].split()[0])
            except ValueError:
                print("Error al leer cdiva.")

        elif line == "!Center":
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    data["centers"].append([float(x) for x in lines[j].strip().split()])
                except ValueError:
                    print(f"Error al leer coordenadas del centro en línea: {lines[j]}")
                j += 1

        elif line == "!particle semiaxis x y z in nm":
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    semiaxis_values = [float(x) for x in lines[j].strip().split()]
                    data["R"].append(semiaxis_values[0]) 
                except ValueError:
                    print(f"Error al leer semiejes en línea: {lines[j]}")
                j += 1

        elif line == "! segment lengths":
            try:
                data["lseg"] = float(lines[i+1].split()[1])
            except ValueError:
                print("Error al leer lseg.")

        elif line == "!PBC PBC xmin xmax ymin ymax zmin zmax, 1=yes, 2=wall, 0=bulk":
            try:
                data["PBC"] = ([float(x) for x in lines[i+1].split()[1:]])
            except ValueError:
                print(f"Error al leer coordenadas PBC")

        i += 1
    return data

def path_carpeta(folder_init,n):
    x = os.path.dirname(folder_init)
    for i in range(n-1):
        x = os.path.dirname(x)
    return x

def extract_R_bin(definitions_path):
    R1_np, R2_np = None, None
    lines = read_DEF(definitions_path)
    size_index = None
    for i, line in enumerate(lines):
        if line.strip() == "!particle semiaxis x y z in nm":
            size_index = i + 1
            break

    if size_index is not None:
        try:
            R1_np = float(lines[size_index].split()[0]) 
            j = 0
            while R2_np == None:
                if float(lines[size_index+ j].split()[0]) != R1_np:
                    R2_np = float(lines[size_index + j].split()[0])
                j += 1

        except ValueError:
            R2_np = R1_np
    
    return R1_np, R2_np

def extract_cdiva(definitions_path):
    cdiva = None
    lines = read_DEF(definitions_path)
    size_index = None
    for i, line in enumerate(lines):
        if line.strip() == "!cdiva":
            size_index = i + 1
            break

    if size_index is not None:
        cdiva = float(lines[size_index]) 

    return cdiva

def extract_references(csv_path):
    """Extrae los datos del archivo referencias.csv."""
    references = []
    with open(csv_path, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            references.append(row)
    return references

def update_particle_sizes(lines, gamma, R_np, n1_k_bin, n2_k_bin):
    size_index = None
    nseg = []
    for i, line in enumerate(lines):
        if line.startswith("!chains lenght"):
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    nseg.append([float(x) for x in lines[j].strip().split()][0])
                except ValueError:
                    print(f"Error al leer semiejes en línea: {lines[j]}")
                j += 1

        elif line.startswith("! segment lengths"):
            size_index = i+1
            lseg = float(lines[size_index].split()[1])

    from scipy.optimize import fsolve

    def cubic_equation(x, a, b):
        return x**3 + a*x**2 - b
    def solve_cubic(a, b, x0=1.0):
        root = fsolve(cubic_equation, x0, args=(a, b))
        return root[0]

    size_index = None
    for i, line in enumerate(lines):
        if line.strip() == "!particle semiaxis x y z in nm":
            size_index = i + 1
            break
    if size_index is not None:
        try:
            for n in np.arange(0,n1_k_bin):
                lines[size_index] = f"{R_np} {R_np} {R_np}\n"
            D1 = 2*R_np
            l = lseg*np.cos(68*np.pi/180/2)
            h1 = (nseg[0]*l+0.2)
            h2 = (nseg[-1]*l+0.2)
            lamda1 = 2*h1/D1
            a = 6*h2/2
            b = np.power(gamma*R_np,3) *(1 + 3*lamda1)

            factor = solve_cubic(a, b)

            new_size = factor.round(2)
            for n in np.arange(0,n2_k_bin):
                lines[size_index + n1_k_bin + n] = f"{new_size} {new_size} {new_size}\n"
        except (ValueError, IndexError):
            print("Find an error in reading or updating the NP size")
    else:
        print("Couldn`t find the NP size seccion in DEFINITIONS lines.")
    
    return lines

def update_R1(DEF, n1_k_bin, R1):
    if os.path.exists('DEFINITIONS_backup.txt'):
        shutil.copy('DEFINITIONS_backup.txt', os.path.join(os.getcwd(), "DEFINITIONS.txt"))
    DEF = "DEFINITIONS.txt"
    lines = read_DEF(DEF)

    size_index = None
    for i, line in enumerate(lines):
        if line.strip() == "!particle semiaxis x y z in nm":
            size_index = i + 1
            break
    for n in np.arange(0,n1_k_bin):
        lines[size_index+n] = f"{R1} {R1} {R1}\n"

    write_DEF(DEF, lines)
    shutil.copy(DEF, os.path.join(os.getcwd(), "DEFINITIONS_backup.txt"))

def generate_references_csv(references, output_folder, delta_value, R, dimx, dimy, dimz, label, name_bin, gamma_value):
    references_path = os.path.join(output_folder, "tot_references.csv")
    existing_dict = defaultdict(set)  # Usar defaultdict para automáticamente manejar los conjuntos de gamma

    # Leer el archivo existente si ya existe
    if os.path.exists(references_path):
        with open(references_path, mode='r', newline='') as file:
            reader = csv.reader(file)
            header = next(reader, None)  # Guardamos la cabecera
            for row in reader:
                row_base = tuple(row[:-1])  # Todo menos gamma
                gamma_field = row[-1].strip()  # Campo gamma
                gamma_list = [float(g) for g in gamma_field.split(',')]
                existing_dict[row_base].update(gamma_list)

    # Procesar nuevas referencias
    for row in references[1:]:
        pos_tuple = tuple(f"{float(val):.10f}" for val in row[1:4])
        hash_input = ",".join(pos_tuple)
        key_value = f"key_{hashlib.md5(hash_input.encode()).hexdigest()[:8]}"

        # Construir la tupla base sin el valor de gamma
        row_base = tuple([label, R] + row + [delta_value, dimx, dimy, dimz, key_value])

        # Agregar el valor de gamma al conjunto en existing_dict, evitando duplicados
        existing_dict[row_base].add(gamma_value)

    # Reescribir el archivo con todos los datos unificados
    with open(references_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Escribir la cabecera
        writer.writerow(['#part', 'radius [nm]'] + references[0] + ['delta', 'dimx', 'dimy', 'dimz', 'key', 'gamma'])
        
        # Escribir las filas con las claves y valores de gamma
        for row_base, gamma_set in existing_dict.items():
            # Asegurarse de que gamma_set contiene floats antes de formatear
            gamma_str = ','.join(f"{g:.3f}" for g in sorted(gamma_set))
            writer.writerow(row_base + (gamma_str,))

def gamma_calc(definitions_path):
    lines = read_DEF(definitions_path)
    size_index = None
    R2_np = None
    nseg_var = []
    nseg = []
    for i, line in enumerate(lines):
        if line.startswith("!chains lenght"):
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    nseg_var.append([float(x) for x in lines[j].strip().split()][0])
                except ValueError:
                    print(f"Error al leer semiejes en línea: {lines[j]}")
                j += 1

        if line.strip() == "!properties of ligand chains":
            nseg = float(lines[i+1].split()[1])

        if line.startswith("! segment lengths"):
            size_index = i+1
            lseg = float(lines[size_index].split()[1])
        if line.strip() == "!particle semiaxis x y z in nm":
            size_index = i + 1
            R1_np = float(lines[size_index].split()[0])
            j = 0
            while R2_np == None:
                try:
                    if R1_np != float(lines[size_index + j].split()[0]):
                        R2_np = float(lines[size_index + j].split()[0])
                except ValueError:
                    R2_np = R1_np
                j += 1
            size_index = None

    if len(nseg_var)>1:
        chain_lenght = {"part1": nseg_var[0], "part2": nseg_var[-1]}
    else: chain_lenght = {"part1": None, "part2": None}

    for label in ['part1','part2']:
        if chain_lenght[label] == None:
            chain_lenght[label] = nseg

    l = lseg*np.cos(68*np.pi/180/2)
    h1 = (chain_lenght['part1']*l+0.2); h2 = (chain_lenght['part2']*l+0.2)
    lamda_fact1 = 2*h1; lamda_fact2 = 2*h2
    D1 = 2*R1_np ;D2 = 2*R2_np
    gamma = R2_np*np.power(1+3*lamda_fact2/D2,(1./3.)) / (R1_np*np.power((1+3*lamda_fact1/D1),(1./3.)))
    return gamma

def join_F_csv(folder, name, bin_true):
    import os, re, glob
    import pandas as pd

    ruta_archivos = folder
    pattern = re.compile(rf"{name}_results_(.+)\.csv$")
    csv_DEFs = [f for f in glob.glob(os.path.join(ruta_archivos, "*.csv")) if pattern.search(os.path.basename(f))]

    data_merged = None

    for i, file in enumerate(csv_DEFs):
        match = pattern.search(os.path.basename(file))
        if not match:
            continue
        f_name = match.group(1)

        df = pd.read_csv(file)

        # Separar primeras 6 columnas + F_value
        first_cols = df.columns[:6].tolist()
        df = df[first_cols + ["F_value"]]
        df.rename(columns={"F_value": f"{f_name}"}, inplace=True)

        if data_merged is None:
            data_merged = df
        else:
            data_merged = pd.merge(data_merged, df, on=first_cols, how='outer')

        os.remove(file)

    # Encabezado artificial
    header_row = [""] * 6 + [f"{key}" for key in data_merged.columns[6:]]
    df_header = pd.DataFrame([header_row], columns=data_merged.columns)

    # Concatenar encabezado y datos
    df_final = pd.concat([df_header, data_merged], ignore_index=True)
    df_final.columns = [col.strip().replace("'", "").replace("(", "").replace(")", "").replace(",", "") for col in df_final.columns]

    data_final = df_final.drop(index=0, errors='ignore')
    output_DEF = os.path.join(ruta_archivos, f"{name}_results_output.csv")
    data_final.to_csv(output_DEF, index=False)

def join_F_csv_ref(folder, name):
    import os, re, glob
    import pandas as pd

    ruta_archivos = folder
    pattern = re.compile(rf"{name}_references_(.+)\.csv$")
    csv_DEFs = [f for f in glob.glob(os.path.join(ruta_archivos, "*.csv")) if pattern.search(os.path.basename(f))]

    data_merged = None

    for i, file in enumerate(csv_DEFs):
        match = pattern.search(os.path.basename(file))
        if not match:
            continue
        f_name = match.group(1)

        df = pd.read_csv(file)

        # Separar primeras 6 columnas + F_reference
        first_cols = df.columns[:6].tolist()
        df = df[first_cols + ["F_reference"]]
        df.rename(columns={"F_reference": f"{f_name}_reference"}, inplace=True)

        if data_merged is None:
            data_merged = df
        else:
            data_merged = pd.merge(data_merged, df, on=first_cols, how='outer')

        os.remove(file)

    header_row = [""] * 6 + [col for col in data_merged.columns[6:]]
    df_header = pd.DataFrame([header_row], columns=data_merged.columns)

    df_final = pd.concat([df_header, data_merged], ignore_index=True)
    df_final.columns = [col.strip().replace("'", "").replace("(", "").replace(")", "").replace(",", "") for col in df_final.columns]

    data_final = df_final.drop(index=0, errors='ignore')
    output_DEF = os.path.join(ruta_archivos, f"{name}_references_output.csv")
    data_final.to_csv(output_DEF, index=False)

def vol_tot_bin(name,R1,R2,l_pol_1,l_pol_2,sigma_1,sigma_2,V_pol):
    pi = math.pi
    Vol_NP_1 = pi*(4./3.)*np.power(R1,3)
    Vol_NP_2 = pi*(4./3.)*np.power(R2,3)
    A_1 = 4*pi*np.power(R1,2)
    A_2 = 4*pi*np.power(R2,2)

    if name == "NaCl":
        N1 = 4
        N2 = 4
    elif name == "CsCl":
        N1 = 1
        N2 = 1
    elif name == "CaCu5":
        N1 = 1
        N2 = 5

    elif name == "AlB2":
        N1 = 1
        N2 = 2

    elif name == "Li3Bi":
        N1 = 0.5
        N2 = 1.5

    elif name == "AuCu":
        N1 = 2
        N2 = 2

    elif name == "Cu3Au":
        N1 = 1
        N2 = 3

    elif name == "NaZn13":
        N1 = 1
        N2 = 12

    elif name == "CaB6":
        N1 = 1
        N2 = 6

    elif name == "Fe4C":
        N1 = 1
        N2 = 4

    elif name == "bccAB6":
        N1 = 2/2**3
        N2 = 12/2**3

    elif name=="MgZn2":
        N1 = 4
        N2 = 8

    v_tot = N1*Vol_NP_1+N1*float(sigma_1)*A_1*V_pol*float(l_pol_1) + N2*Vol_NP_2+N2*float(sigma_2)*A_2*V_pol*float(l_pol_2)

    return v_tot

