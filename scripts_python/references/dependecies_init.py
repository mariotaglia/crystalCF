import numpy as np
import os
import shutil

def read_DEF(file_path):
	"""Extract the lines from DEF."""
	with open(file_path, 'r') as f:
		lines = f.readlines()
	return lines

def write_DEF(file_path, lines):
	"""Export DEF with lines."""
	with open(file_path, 'w') as f:
		f.writelines(lines)

def extract_definitions(definitions_path):
    data = {
        "centers": [],
        "R": [],
    }
    lines = read_DEF(definitions_path)
    for i, line in enumerate(lines):
        line = lines[i].strip()
        if line == "!Center":
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

        i += 1
    return data

def edit_def_bin(DEF, n1, n2, k_bin, name):
	lines = read_DEF(DEF)

	sections_info = [
		("!cdiva",1,1),
		("!gama",1,1),
		("!PBC PBC xmin xmax ymin ymax zmin zmax, 1=yes, 2=wall, 0=bulk",1,1),
		("! number of particles", 1, 1),
		("!Center", 1, k_bin*(n1+n2)),
		("!particle semiaxis x y z in nm", 1, k_bin*(n1+n2)),
		("! Rotation", 3, k_bin*(n1+n2)),
		("! coverage", 1, k_bin*(n1+n2)),
		("!Surface-polymer atraction", 1, k_bin*(n1+n2))
	]
	if name=="MgZn2":
		modified_lines = lines.copy()
		particles_to_keep = [2,3,6,7,8,9,10,11,12]
		num_particles = len(particles_to_keep)

		sections_found = []
		for key, lines_per_particle, tot_particles in sections_info:
			for i, line in enumerate(lines):
				if line.strip() == key:
					start_index = i + 1
					sections_found.append((key, start_index, lines_per_particle, tot_particles))
					break

		for key, start_index, lines_per_particle, tot_particles in sorted(sections_found, key=lambda x: x[1], reverse=True):
			if key == "! number of particles":
				modified_lines[start_index] = f"{num_particles}\n"
			elif key == "!cdiva":
				modified_lines[start_index] = f"{str((8/3)**0.5/2)}\n"
			elif key == "!gama":
				continue
			elif key == "!PBC PBC xmin xmax ymin ymax zmin zmax, 1=yes, 2=wall, 0=bulk":
				modified_lines[start_index] = f"PBC 1 1 1 1 3 3\n"
			else:
				new_block = []
				for part in particles_to_keep:
					actual_index = start_index + (part - 1) * lines_per_particle
					new_block.extend(modified_lines[actual_index:actual_index + lines_per_particle])

				modified_lines[start_index:start_index + tot_particles * lines_per_particle] = []
				modified_lines[start_index:start_index] = new_block

		output_DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
		with open(output_DEF, "w") as f:
			f.writelines(modified_lines)

		data = extract_definitions("DEFINITIONS.txt")
		centers = data.get("centers", [])
		lines = read_DEF('DEFINITIONS.txt')

		sections_info = {
		"!Center": num_particles,  
		}

		modified_lines = lines.copy()
		for i, line in enumerate(modified_lines):
			if line.strip() == "!Center":
				start_index = i + 1
				break

		for j in range(num_particles):
			pos1, pos2, pos3 = centers[j]
			modified_lines[start_index+j] = f"{pos1} {pos2} {(pos3-1/4)*2}\n"

		output_DEF = os.path.join(os.getcwd(), "DEFINITIONS.txt")
		with open(output_DEF, "w") as f:
			f.writelines(modified_lines) 

def list_reflexion(name):
	list_reflexion = [
	("NaCl",[1,1]),
	("CsCl",[1,1]),
	("MgZn2",[2,7])
	]

	for name_bin, n_k_list in list_reflexion:
		if name_bin == name:
			list_reflexion = n_k_list

	return list_reflexion
