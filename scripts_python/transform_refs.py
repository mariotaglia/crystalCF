import subprocess
import csv

def read_DEF(file_path):
    """Extract the lines from DEF."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return lines

def x_center(center, Nref, delta, dimx, aL):
    x = center * delta * dimx
    x_grid = x % delta
    x0 = delta * (Nref - 1) / 2
    x_ref = x0 + x_grid
    frac_ref_x = x_ref / aL
    return frac_ref_x

def z_center(center, Nref, delta, dimz, cdiva, aL):
    z = center * delta * cdiva * dimz
    z_grid = z % (delta * cdiva)
    z0 = delta * cdiva * (Nref - 1) / 2
    z_ref = z0 + z_grid
    frac_ref_z = z_ref / (aL * cdiva)
    return frac_ref_z

def calculate_center_ref(Nref, centers, dimx, dimy, dimz, delta, cdiva, PBC):
    aL = Nref * delta
    center_ref_list = []
    for center in centers:
        if PBC[0] == 3 and PBC[1] == 3:
            if center[0] == 0 or center[0] == 1:
                frac_ref_x = center[0]
            else:
                frac_ref_x = x_center(center[0], Nref, delta, dimx, aL)
        else:
            frac_ref_x = x_center(center[0], Nref, delta, dimx, aL)

        if PBC[2] == 3 and PBC[3] == 3:
            if center[1] == 0 or center[1] == 1:
                frac_ref_y = center[1]
            else:
                frac_ref_y = x_center(center[1], Nref, delta, dimy, aL)
        else:
            frac_ref_y = x_center(center[1], Nref, delta, dimy, aL)
        
        if PBC[4] == 3 and PBC[5] == 3:
            if center[2] == 0 or center[2] == 1:
                frac_ref_z = center[2]
            else:
                frac_ref_z = z_center(center[2], Nref, delta, dimz, cdiva, aL)
        else:
            frac_ref_z = z_center(center[2], Nref, delta, dimz, cdiva, aL)
        
        center_ref_list.append([frac_ref_x, frac_ref_y, frac_ref_z])
    
    return center_ref_list

def process_positions(center_ref_list):
    pos_out = []
    counts = []
    tol = 1e-4
    
    for center_ref in center_ref_list:
        found = False
        for i, existing_ref in enumerate(pos_out):
            if all(abs(a - b) < tol for a, b in zip(center_ref, existing_ref)):
                counts[i] += 1
                found = True
                break
        if not found:
            pos_out.append(center_ref)
            counts.append(1)
    
    with open('references.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        
        writer.writerow(["Count", "Pos 1", "Pos 2", "Pos 3"])
        
        for count, pos in zip(counts, pos_out):
            writer.writerow([count] + pos)

    return pos_out, counts

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
        "nseg": None
    }
    lines = read_DEF(definitions_path)
    for i, line in enumerate(lines):
        line = lines[i].strip()
        parts = line.split()
        if len(parts) == 2 and parts[0] in data:
            try:
                data[parts[0]] = float(parts[1])
            except ValueError:
                print(f"Error while reading {parts[0]}: {parts[1]}.")

        elif line == "! number of particles":
            try:
                data["num_particles"] = float(lines[i+1].split()[0])
            except ValueError:
                print("Error while reading the number of particles.")
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
                    print(f"Error while reading the !center: {lines[j]}")
                j += 1

        elif line == "!particle semiaxis x y z in nm":
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    semiaxis_values = [float(x) for x in lines[j].strip().split()]
                    data["R"].append(semiaxis_values[0]) 
                except ValueError:
                    print(f"Error while reading particle semiaxis: {lines[j]}")
                j += 1
        elif line == "!properties of ligand chains":
            try:
                data["nseg"] = float(lines[i+1].split()[1])
            except ValueError:
                print("Error while reading \'long\'")
        elif line == "! segment lengths":
            try:
                data["lseg"] = float(lines[i+1].split()[1])
            except ValueError:
                print("Error while reading \'lseg\'")
        i += 1
    return data

def main():
    definitions_path = "DEFINITIONS.txt"
    data = extract_definitions(definitions_path)

    dimx = int(data.get("dimx"))
    dimy = int(data.get("dimy"))
    dimz = int(data.get("dimz"))
    delta = float(data.get("delta"))
    cdiva = float(data.get("cdiva"))
    centers = data.get("centers", [])

    center_ref_list = calculate_center_ref(centers, dimx, dimy, dimz, delta, cdiva)
    pos_out, counts = process_positions(center_ref_list)

    with open('references.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        
        writer.writerow(["Count", "Pos 1", "Pos 2", "Pos 3"])
        
        for count, pos in zip(counts, pos_out):
            writer.writerow([count] + pos)

    print("pos_out:")
    for pos in pos_out:
        print(f"{pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}")

if __name__ == '__main__':
    main()
