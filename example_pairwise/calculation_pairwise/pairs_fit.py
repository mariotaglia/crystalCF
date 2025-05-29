import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os 
from scipy import interpolate

fig,ax = plt.subplots(2,figsize=(11,10))

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

ydata = np.zeros(0)
xdata = np.zeros(0)

Radios = []
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
         Dsopm[:] = Ds[:] - opm(Rs[0],12) - opm(rB,12) - 4.0 # 4.0 is the shift in MOLT calculation
         Fnorm = -Fs/F0(Rs[0],rB) 

         ax[1].scatter(Dsopm,Fnorm,s=3)
         ax[0].plot(Ds,Fs)
  
         Radios.append("R: "+str(Rs[0]) + " nm")
     
         ydata = np.append(ydata,Fnorm)
         xdata = np.append(xdata,Dsopm)

xydata = np.transpose(np.stack((xdata,ydata)))

ax[0].legend(Radios,loc='lower right',fontsize='x-small')         

xydata = xydata[xydata[:, 0].argsort()]

tck = interpolate.splrep(xydata[:,0], xydata[:,1], k=2, s=5)



xnew = np.linspace(np.min(xdata), 13, 50)
ynew = interpolate.splev(xnew, tck, der=0)


xynew = np.transpose(np.stack((xnew,ynew)))


#ax.scatter(xydata[:,0],xydata[:,1])
ax[1].plot(xnew, ynew, label = 'Fit')

ax[0].set_xlabel('D / nm')
ax[1].set_xlabel('D - Dopm / nm')
ax[0].set_ylabel('$F_{pair} / k_BT$') 
ax[1].set_ylabel('$F_{pair}/F_0$')

np.savetxt('fitpairL12.dat',xynew)
plt.show()


