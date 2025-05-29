import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os 

fig,ax = plt.subplots()

def opm(r,L):
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

rB = 2
for filename in os.listdir("."):
    if filename.endswith(".dat") and filename.startswith("pair"):
 
         data = pd.read_csv(filename, sep=" ", header=0, names=["R","D","case", "F"])
 
         Ds = data["D"].to_numpy()
         Fs = data["F"].to_numpy()
         Rs = data["R"].to_numpy()

         Fs[:] = Fs[:] - Fs[-1]

         print(Rs[0],np.min(Fs))

         Dsopm = np.zeros(np.size(Ds))
         Fnorm = np.zeros(np.size(Ds))
         Dsopm[:] = Ds[:] - opm(Rs[0],12)-opm(rB,12)
         Fnorm = -Fs/F0(Rs[0],rB) 

         ax.scatter(Dsopm,Fnorm)
#         ax.plot(Ds,Fs)


plt.show()

