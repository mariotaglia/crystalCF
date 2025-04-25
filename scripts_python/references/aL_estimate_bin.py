import numpy as np
import math
import sys
import os

if len(sys.argv) < 4:
    raise ValueError("You need to set name_bin, R1 and R2 params")

name = str(sys.argv[1])
R1 = float(sys.argv[2])
R2 = float(sys.argv[3])

sigma = 5.85 #1/nm^2
pi = math.pi
Vol_NP_1 = pi*(4/3)*R1**3
Vol_NP_2 = pi*(4/3)*R2**3
A_1 = 4*pi*R1**2
A_2 = 4*pi*R2**2
l_pol = 12

if name == "NaCl":
	V_pol = 0.0285 #nm^3
	N1 = 2
	N2 = 2
elif name == "CsCl":
	V_pol = 0.028 #nm^3
	N1 = 1
	N2 = 1
elif name == "MgZn2":
	V_pol = 0.013
	N1 = 4
	N2 = 8
elif name == "CaCu5":
	V_pol = 0.031
	N1 = 1
	N2 = 5

elif name == "AlB2":
	V_pol = 0.035
	N1 = 1
	N2 = 2

elif name == "Li3Bi":
	V_pol = 0.071
	N1 = 2
	N2 = 6

elif name == "AuCu":
	if R2<=R1*0.88/2:
		V_pol = 0.030
		N1 = 2
		N2 = 2
	else:
		V_pol = 0.035
		N1 = 2
		N2 = 2

elif name == "Cu3Au":
	V_pol = 0.038
	N1 = 1
	N2 = 3

elif name == "NaZn13":
	V_pol = 0.033
	N1 = 1
	N2 = 12

elif name == "CaB6":
	V_pol = 0.031
	N1 = 1
	N2 = 6

if name=="MgZn2":
	if R2<=R1*1.12/2:
		R = R1*0.88/2
		V_pol = 0.014
		A_2 = 4*pi*R**2
		Vol_NP_2 = pi*(4./3.)*R**3

	elif R2<=R1*1.48/2:
		V_pol = 0.0115

	else:
		V_pol = 0.010

	aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+ N2*Vol_NP_2+N2*sigma*A_2*V_pol*l_pol)**(1./3.)

elif name == "NaZn13":
	aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+ N2*Vol_NP_2+N2*sigma*A_2*V_pol*l_pol)**(1./3.) *2

else:
	aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+ N2*Vol_NP_2+N2*sigma*A_2*V_pol*l_pol)**(1./3.)

def output():
	print(str(aL_min))
if __name__ == '__main__':
	output()
