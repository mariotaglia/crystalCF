import numpy as np
import os.path
import sys

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



def read_radii():
    count = 0
    NNN = 0

    with open("DEFINITIONS.txt") as fp:
        l = fp.readlines()
        for line in l:
            count += 1
            if 'number of particles' in line:
                pos = count
        NNN = int(l[pos])
 
        radii = np.zeros(NNN)
        pos = pos + 2 + NNN
        for i in range(NNN):
            line = l[pos+i+1]
            radii[i]=line.split()[0]
    return NNN, radii

######################################################
# MODIFY DEFINITIONS.txt
######################################################

def cdiva_calc(filename,gamma):
    name = filename.split(".")[0]
    flag = 0
    if name == "CaCu5":
        if gamma < 2.0 / np.sqrt(3.0) - 1:
            a_fac = 1.0
            c_fac = 1.0
        elif (gamma >= 2.0 / np.sqrt(3.0) - 1) and (gamma < (1 + 2 * np.sqrt(19)) / 15.0):
            a_fac = np.sqrt(3) * (1 + gamma) * 0.5
            c_fac = 1.0
        elif (gamma >= (1 + 2 * np.sqrt(19)) / 15.0) and (gamma<1 / (4.0/np.sqrt(3)-1.0)):
            a_fac = np.sqrt(3) * (1 + gamma)*0.5
            c_fac = 0.5 * np.sqrt(15.0*gamma**2 - 2*gamma - 1)
        else:
            a_fac = 2*gamma
            c_fac = np.sqrt(8.0/3.0)*gamma
        flag = 1

    if name == "AlB2":
        if gamma < np.sqrt(7.0/3.0)-1:
            a_fac = 1.0
            c_fac = 1.0
        elif gamma > 2.0/3.0:
            a_fac = np.sqrt(3.0)*gamma
            c_fac = 1.0
        elif (gamma > 1/np.sqrt(3.0)) and (gamma <= 2.0/3.0):
            a_fac = np.sqrt(3.0)*gamma
            c_fac = np.sqrt(-3.0*gamma**2+2.0*gamma+1.0)
        else:
            c_fac = np.sqrt(gamma**2+2*gamma-1.0/3.0)
            a_fac = 1.0
        flag = 1

    if name == "AuCu":
        if (0 < gamma) and (gamma < np.sqrt(3.0)-1):
            a_fac = np.sqrt(2)
            c_fac = 1

        elif (gamma >= np.sqrt(3.0)-1) and (gamma <= 1.0):
            a_fac = np.sqrt(2)
            c_fac = np.sqrt(2)*np.sqrt(0.5*gamma**2+gamma-0.5)

        else:
            raise ValueError('Parameter gamma is outside allowed range')
        flag = 1

    if name == "Cu3Au":
        if gamma < np.sqrt(2.0)-1:
            c_fac = 1
        else:
            c_fac = (1+gamma)/np.sqrt(2.0)
        a_fac = c_fac
        flag = 1

    if name == "MgZn2":
        a_fac = 1
        c_fac = 1.632993162
        flag = 1

    if name == "FCC":
        a_fac = 1
        c_fac = np.sqrt(2.0)
        flag = 1

    if name == "BCC":
        a_fac = 1
        c_fac = 1
        flag = 1

    if name == "CsCl":
        a_fac = 1
        c_fac = 1
        flag = 1
 
    if name == "Li3Bi":
        a_fac = 1
        c_fac = 1
        flag = 1
        
    if flag == 0:
        print("Structure unknown")
        exit()

    cdiva = c_fac/a_fac
    if cdiva >= 2 or cdiva <= 0.5:
        print(f"warning, cdiva<=0.5 or cdiva>=2.0 for gamma {gamma:.2f}")
    return cdiva


def modify_definitions(filename, rA, rB, real_gamma):

    cdiva = cdiva_calc(filename,real_gamma)
    
    with open(filename, "r") as f:
        lines = f.readlines()
    with open("DEFINITIONS.txt", "w") as f:
        for line in lines:
            if ("dumpcluster" not in line) and  ("cutoffcluster" not in line) and ("cluster_same" not in line):
                line = line.replace("_rA",str(rA))
                line = line.replace("_rB",str(rB))
                line = line.replace("_cdiva",str(cdiva))
                f.write(line)

        f.write("dumpcluster 2 \n")
        f.write("cutoffcluster 15 \n")
        f.write("cluster_same 0 \n")


def modify_definitions2(aL):
    with open("DEFINITIONS.txt", "r") as f:
        lines = f.readlines()
    with open("DEFINITIONS.txt", "w") as f:
        for line in lines:
            line = line.replace("_aL",str(aL))
            f.write(line)

####################################################

def calc_F(radii):

    os.system("~/develop/crystalCF/crystalCF > output.dat")
    flag = 0
    
#####################################################
# Read distances
######################################################

    if os.path.isfile('distances.dat') == False:
        print('Failed to find distances.dat, increase packmin')
        exit()

    with open("distances.dat") as fp:
        l = fp.readlines()
        numdists = int(l[0])
        if numdists == 0: # particles too far away
            return 0.0

#        print("Number of pairs to consider:", numdists)
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
            
#            exit()
        if (dists[i] - rA - rB < 0): 
            print("Cores in contact") 
            FAB = 1e100
#            exit()
        Fpair += -FAB*FAB0/2.*weights[i]   

#        print("Pairwise additive free-energy:", Fpair)

#    with open("F_pairwise.dat","w") as f:
#        f.write(str(1)+"   "+str(Fpair)+"\n")

    return Fpair

#####################################################

print()
print('------------------------------------------------------------------------------------------------------------')
print('Calculation of F vs gamma from additive free energy + packing')
print('Needs ~/develop/crystalCF/crystalCF"')
print('Need fitpairL12.dat calculated with pairs_fit.py')
print('Command line arguments:')
print('file.txt : file to use as DEFINITIONS.txt')
print('packmin minimum packing fraction to use in calculation')
print('------------------------------------------------------------------------------------------------------------')
print()

###########################################################
# READ number of particles and their sizes from DEFINITIONS
###########################################################

filename=sys.argv[1]
packmin=float(sys.argv[2])

if os.path.isfile(filename) == False:
    print('Failed to find' + filename)
    exit()

folder = filename.split(".")[0]


os.system("mkdir "+folder)
os.chdir(folder)

##########################################################
# Loop over gamma
##########################################################

rB = 2. # size of the largest particle
sigma = 5.85 # surface coverage in nm^-2
vpol = 0.028 # bead (methylene) volume in nm^3
longc = 12 # ligand chain length 
pi = 3.14159
jj = 0

gamma_core_list = np.linspace(.2,1.0,17) # value of gamma core to use
F_list = np.zeros(gamma_core_list.size)
pack_list = np.zeros(gamma_core_list.size)
real_gamma_list = np.zeros(gamma_core_list.size)

for gamma_core in gamma_core_list:


    rA = rB * gamma_core # size of smallest particle
    real_gamma = ropm(rA,12)/ropm(rB,12)

    folder = 'gamma_'+str(gamma_core)
    os.system("mkdir "+folder)
    os.chdir(folder)
    os.system("cp ../../"+filename+" .")
    os.system("cp ../../fitpairL12.dat .")

    modify_definitions(filename, rA, rB, real_gamma) # changes radii in filename and renames to DEFINITIONS.txt

    NNN, radii = read_radii() # reads radii list and number of particles

    packing_list = np.linspace(packmin,1.0,100)

    FPACK = np.ones(packing_list.size)*1e100

    ii = 0 # counter

#    print(ropm(rA,12)/ropm(rB,12), real_gamma_list[ii], ii)
    for pack in packing_list: # loop over values of packing
        print(pack)

        folder = 'pack_'+str(pack)
        os.system("mkdir "+folder)
        os.chdir(folder)
        os.system("cp ../DEFINITIONS.txt .")
        os.system("cp ../fitpairL12.dat .")


#### Determine Packing #####
        particle_vol = 0.0
        particle_surf = 0.0
        for i in range(NNN):
            particle_vol += 4./3.*pi*np.power(radii[i],3)
            particle_surf += 4.*pi*np.power(radii[i],2)
        total_vol = particle_vol + particle_surf*sigma*longc*vpol # total particle volume
#############################


        aL = np.power(total_vol/pack, 1./3.) # lattice parameter
        modify_definitions2(aL) # changes radii in filename and renames to DEFINITIONS.txt
        FPACK[ii] = calc_F(radii)
        ii = ii+1

        os.chdir('..') # go back to gamma_

        
    F_list[jj] = np.min(FPACK)
    real_gamma_list[jj] = real_gamma
    pack_list[jj] = packing_list[np.argmin(FPACK)]
    if pack_list[jj] == packmin:
        print("Decrease packmin")
        exit()
    print(gamma_core_list[jj])
    jj = jj + 1
    os.chdir('..') # go back to BCC



for i in range(pack_list.size):
    print(gamma_core_list[i], real_gamma_list[i],  pack_list[i], F_list[i]/NNN)


data = np.column_stack((real_gamma_list, pack_list, F_list/NNN))

# Save to CSV
np.savetxt('output.csv', data, delimiter=',', header='Gamma,Packing,F', comments='')

