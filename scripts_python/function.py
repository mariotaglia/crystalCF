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
        "cell part": [], "list delta part": {}, "num cell part": {},
        "num cell bin": None, "aL cell bin factor": None, 
        "aL cell part factor": {}, "flag generate energy vs aL curves": None,
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
        elif line == "!aL cell bin factor":
            data["aL cell bin factor"] = eval(lines[i + 1].strip('\n'))
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
        elif line == "!aL cell part factor":
            j = i + 1
            while j < len(lines) and not lines[j].startswith("!"):
                parts = lines[j].split(maxsplit=1)
                if len(parts) > 1 and parts[0] in data["cell part"]: 
                    key = parts[0]
                    if key in data["cell part"]:
                        data["aL cell part factor"][key] = eval(parts[1].strip('\n'))
                j += 1
            i = j - 1
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
        "nseg": None,
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
        elif line == "!properties of ligand chains":
            try:
                data["nseg"] = float(lines[i+1].split()[1])
            except ValueError:
                print("Error al leer nseg.")
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
    for i, line in enumerate(lines):
        if line.strip() == "!properties of ligand chains":
            size_index = i+1
            n_seg = float(lines[size_index].split()[1])
            size_index = None
        if line.startswith("! segment lengths"):
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
            D = 2*R_np
            lamda = 2*n_seg*lseg/D
            a = 6*n_seg*lseg/2
            b = (gamma*R_np)**3 *(1 + 3*lamda)

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

def generate_references_csv(references, output_folder, delta_value, R, dimx, dimy, dimz, label, name_bin, k_aL):
    references_path = os.path.join(output_folder, "tot_references.csv")
    existing_rows = []  # Lista para mantener el orden y evitar duplicados

    # Leer el archivo existente si ya existe
    if os.path.exists(references_path):
        with open(references_path, mode='r', newline='') as file:
            reader = csv.reader(file)
            header = next(reader, None)  # Guardamos la cabecera
            existing_rows = [tuple(row) for row in reader]  # Convertir a tuplas para comparación exacta

    new_rows = []

    for row in references[1:]:
        pos_tuple = tuple(f"{float(val):.10f}" for val in row[1:4])
        hash_input = ",".join(pos_tuple)
        key_value = f"key_{hashlib.md5(hash_input.encode()).hexdigest()[:8]}"
        new_row = tuple([label, R] + row + [delta_value, dimx, dimy, dimz , key_value])
        if new_row not in existing_rows:  # Solo agregar si la fila no existe
            new_rows.append(new_row)
            existing_rows.append(new_row)  # Agregarla a la lista de comparación

    if new_rows:
        file_exists = os.path.exists(references_path)

        with open(references_path, mode='a', newline='') as file:
            writer = csv.writer(file)

            if not file_exists:  # Si el archivo no existía, escribir la cabecera
                writer.writerow(['#part', 'radius [nm]'] + references[0] + ['delta', 'dimx', 'dimy', 'dimz', 'key'])

            writer.writerows(new_rows)  # Escribir solo las nuevas filas

    return references_path


def gamma_calc(definitions_path):
    lines = read_DEF(definitions_path)
    size_index = None
    R2_np = None
    for i, line in enumerate(lines):
        if line.strip() == "!properties of ligand chains":
            size_index = i+1
            n_seg = float(lines[size_index].split()[1])
            size_index = None
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

    lamda_fact = 2*lseg*n_seg
    D1 = 2*R1_np ;D2 = 2*R2_np
    gamma = R2_np*(1+3*lamda_fact/D2)**(1./3.) / (R1_np*(1+3*lamda_fact/D1)**(1./3.))
    return gamma

def join_F_csv(folder, name, bin_true):
    ruta_archivos = folder
    pattern = re.compile(rf"{name}_results_(.+)\.csv$")

    # Filtrar archivos que cumplen con el patrón
    csv_DEFs = [f for f in glob.glob(os.path.join(ruta_archivos, "*.csv")) if pattern.search(os.path.basename(f))]

    # Diccionario para almacenar los DataFrames y sus respectivos f_name
    data_dict = {}

    for i, file in enumerate(csv_DEFs):
        # Extraer f_name} desde el nombre del archivo
        match = pattern.search(os.path.basename(file))
        if not match:
            continue
        f_name = match.groups()
        # Cargar el CSV
        if i==0:
            df = pd.read_csv(file)
        else:
            df = pd.read_csv(file,usecols=["F_value"])

        df.rename(columns={"F_value": f"{f_name}"}, inplace=True)
        data_dict[f_name] = df

        os.remove(file)

    # Asegurar que todas las columnas de los DataFrames sean las mismas, llenando NaN si falta alguna
    df_combined = pd.concat(data_dict.values(), axis=1, join='outer')
    # Crear la fila extra con los nombres de f_name en la posición correcta
    if bin_true == True:
        header_row = ["", "","", "", "", ""] + [f"{f_name}" for f_name in data_dict.keys()]  # El encabezado con los nombres f_name
    else:
        header_row = ["", "", "","", "", ""] + [f"{f_name}" for f_name in data_dict.keys()]  # El encabezado con los nombres f_name
    
    df_header = pd.DataFrame([header_row], columns=df_combined.columns)

    # Concatenar la fila extra antes del DataFrame combinado
    df_final = pd.concat([df_header, df_combined], ignore_index=True).drop(0, axis=0, errors='ignore')
    df_final.columns = [col.strip().replace("'", "").replace("(", "").replace(")", "").replace(",", "") for col in df_final]
    # Guardar el CSV final
    output_DEF = os.path.join(ruta_archivos, f"{name}_results_output.csv")
    df_final.to_csv(output_DEF, index=False)

def join_F_csv_ref(folder, name):
    ruta_archivos = folder
    pattern = re.compile(rf"{name}_references_(.+)\.csv$")

    # Filtrar archivos que cumplen con el patrón
    csv_DEFs = [f for f in glob.glob(os.path.join(ruta_archivos, "*.csv")) if pattern.search(os.path.basename(f))]

    # Diccionario para almacenar los DataFrames y sus respectivos f_name
    data_dict = {}

    for i, file in enumerate(csv_DEFs):
        # Extraer f_name} desde el nombre del archivo
        match = pattern.search(os.path.basename(file))
        if not match:
            continue
        f_name = match.groups()
        # Cargar el CSV
        if i==0:
            df = pd.read_csv(file)
        else:
            df = pd.read_csv(file,usecols=["F_reference"])

        df.rename(columns={"F_reference": f"{f_name}_reference"}, inplace=True)
        data_dict[f_name] = df

        os.remove(file)

    # Asegurar que todas las columnas de los DataFrames sean las mismas, llenando NaN si falta alguna
    df_combined = pd.concat(data_dict.values(), axis=1, join='outer')
    header_row = ["", "","", "", "", ""] + [f"{f_name}" for f_name in data_dict.keys()]  # El encabezado con los nombres f_name
    df_header = pd.DataFrame([header_row], columns=df_combined.columns)

    # Concatenar la fila extra antes del DataFrame combinado
    df_final = pd.concat([df_header, df_combined], ignore_index=True).drop(0, axis=0, errors='ignore')
    df_final.columns = [col.strip().replace("'", "").replace("(", "").replace(")", "").replace(",", "") for col in df_final]
    # Guardar el CSV final
    output_DEF = os.path.join(ruta_archivos, f"{name}_references_output.csv")
    df_final.to_csv(output_DEF, index=False)

