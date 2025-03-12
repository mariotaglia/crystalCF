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
	V_pol = 0.02325 #nm^3
	N1 = 2
	N2 = 2
elif name == "CsCl":
	V_pol = 0.028 #nm^3
	N1 = 1
	N2 = 1
elif name == "MgZn2":
	V_pol = 0.003
	N1 = 4
	N2 = 8

if name == "NaCl" or "CsCl":
	aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+ N2*Vol_NP_1+N2*sigma*A_2*V_pol*l_pol)**(1./3.)
elif name == "MgZn2":
	aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+ N2*Vol_NP_1+N2*sigma*A_2*V_pol*l_pol)**(1./3.) 


if __name__ == '__main__':
	print(str(aL_min))
