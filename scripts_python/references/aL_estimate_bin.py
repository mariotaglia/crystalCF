import numpy as np
import math
import sys
import os

if len(sys.argv) < 4:
    raise ValueError("You need to set name_bin, R1 and R2 params")

name = str(sys.argv[1])
R = float(sys.argv[2])
sigma = float(sys.argv[3])
l_pol = float(sys.argv[4])
l_extended = 0.128*l_pol+0.2
pi = math.pi
Vol_NP = pi*(4/3)*R**3
A = 4*pi*R**2

if name=="C14":
	N = 12
	V_pol = 0.06
	cdiva = np.sqrt(8/3)/2
	aL_min = (N*Vol_NP+N*sigma*A*V_pol*l_extended)**(1./3.)/cdiva**(1./3.)

def output():
	print(str(aL_min))
if __name__ == '__main__':
	output()