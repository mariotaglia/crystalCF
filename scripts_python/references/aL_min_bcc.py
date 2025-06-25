import numpy as np
import math
import sys
import os

if len(sys.argv) < 4:
    raise ValueError("Debes proporcionar R_np como argumento")

R = float(sys.argv[1])  # Convertir argumento a float
l_pol = float(sys.argv[2])  # Convertir argumento a float
sigma = float(sys.argv[3])

sigma = 5.85 #1/nm^2
pi = math.pi
Vol_NP = pi*(4/3)*R**3
A_NP = 4*pi*R**2
V_pol = 0.03075 #nm^3
N = 2

aL_min = (N*Vol_NP+N*sigma*A_NP*V_pol*l_pol)**(1./3.)

def output():
	print(str(aL_min))
if __name__ == '__main__':
	output()
