import numpy as np
import math
import sys
import os

if len(sys.argv) < 3:
    raise ValueError("Debes proporcionar R1_np y R2_np como argumento")

R1 = float(sys.argv[1])  # Convertir argumento a float
R2 = float(sys.argv[2])  # Convertir argumento a float

sigma = 5.85 #1/nm^2
pi = math.pi
Vol_NP_1 = pi*(4/3)*R1**3
Vol_NP_2 = pi*(4/3)*R2**3
A_1 = 4*pi*R1**2
A_2 = 4*pi*R2**2

V_pol = 0.028 #nm^3
l_pol = 12
N1 = 1
N2 = 1
aL_min = (N1*Vol_NP_1+N1*sigma*A_1*V_pol*l_pol+N2*Vol_NP_1+N2*sigma*A_2*V_pol*l_pol)**(1./3.)

def output():
	print(str(aL_min))
if __name__ == '__main__':
	output()
