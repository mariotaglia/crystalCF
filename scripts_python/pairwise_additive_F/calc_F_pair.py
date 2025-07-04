import numpy as np
import os.path

def ropm(r,L):
    a = 0.09
    b = 0.2
    lb = a*L+b
    ropm = r*np.power((1+3*lb/r),1/3)
    return ropm

def F0(r1,r2):
    a = -211.85
    b = 9.6776
    alfa = -0.55
    reff = np.power((np.power(r1,alfa)+np.power(r2,alfa)),(1/alfa))
    return a*reff+b


print()
print('------------------------------------------------------------------------------------------------------------')
print('Calculation of pairwise additive free energy')
print('Needs ~/develop/crystalCF/crystalCF"')
print('Need fitpairL12.txt calculated with pairs_fit.py')
print('------------------------------------------------------------------------------------------------------------')
print()

######################################################
# READ number of particles and their sizes from DEFINITIONS
######################################################


if os.path.isfile('DEFINITIONS.txt') == False:
    print('Failed to find DEFINITIONS.txt')
    exit()

flag = 0
NNN = 0
count = 0
radii = []

with open("DEFINITIONS.txt") as fp:
    lines = fp.readlines()
    for i, line in enumerate(lines):
        line = lines[i].strip()
        if line == "! number of particles":
            try:
                NNN = int(lines[i+1].split()[0])
            except ValueError:
                print("Error al leer el numero de particulas.")

        elif line == "!particle semiaxis x y z in nm":
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("!"):
                try:
                    semiaxis_values = [float(x) for x in lines[j].strip().split()]
                    radii.append(semiaxis_values[0]) 
                except ValueError:
                    print(f"Error al leer semiejes en línea: {lines[j]}")
                j += 1
            break

######################################################
# MODIFY DEFINITIONS.txt
######################################################


with open("DEFINITIONS.txt", "r") as f:
    lines = f.readlines()
with open("DEFINITIONS.txt_clt", "w") as f:
    for line in lines:
        if ("dumpcluster" not in line) and  ("cutoffcluster" not in line) and ("cluster_same" not in line):
            f.write(line)
    f.write("dumpcluster 2 \n")
    f.write("cutoffcluster 15.0 \n")
    f.write("cluster_same 0 \n")

os.system("mv DEFINITIONS.txt DEFINITIONS.txt_tmp")
os.system("mv DEFINITIONS.txt_clt DEFINITIONS.txt")
os.system("~/develop/crystalCF/crystalCF")
os.system("mv DEFINITIONS.txt_tmp DEFINITIONS.txt")


#####################################################
# Read distances
######################################################

if os.path.isfile('distances.dat') == False:
    print('Failed to find distances.dat')
    exit()

with open("distances.dat") as fp:
    l = fp.readlines()
    numdists = int(l[0])
    #print("Number of pairs to consider:", numdists)
    dists=np.zeros(numdists)
    weights=np.zeros(numdists)
    partA=np.zeros(numdists)
    partB=np.zeros(numdists)
    count = 1
    for i in range(numdists):
        dists[i] = l[count]
        count += 1
        weights[i] = l[count]
        count += 1
        partA[i] = l[count]
        count += 1
        partB[i] = l[count]
        count += 1

######################################################
# Read normalized pair potential
######################################################

if os.path.isfile('fitpairL12.dat') == False:
    print('Failed to find fitpairL12.dat')
    exit()


datanorm = np.loadtxt('fitpairL12.dat')

Dnorm = np.zeros(np.size(datanorm,0))
Fnorm = np.zeros(np.size(datanorm,0))

Dnorm = datanorm[:,0]
Fnorm = datanorm[:,1]

#######################################################
# Calculate pair wise additive F
#######################################################

Fpair = 0.0

for i in range(numdists):
    rA = radii[int(partA[i])-1]
    rB = radii[int(partB[i])-1]
    rA_opm = ropm(rA,12)
    rB_opm = ropm(rB,12)
    FAB0 = F0(rA,rB)

    DAB_norm = dists[i] - rA_opm - rB_opm 

    FAB = np.interp(DAB_norm,Dnorm,Fnorm, left=1e100, right=0.0) # use a linear interpolation to prevent overshooting when the function is close to zero
    if(FAB == 1e100):
        print("Particles too close")
        exit()
    Fpair += -FAB*FAB0/2.*weights[i]   

print("Pairwise additive free-energy:", Fpair)

with open("F_pairwise.dat","w") as f:
    f.write(str(1)+"   "+str(Fpair)+"\n")

