import os
import numpy as np
import matplotlib.pyplot as plt
import sys

vsol = float(sys.argv[1])
temp = float(sys.argv[2])

print("Volume bead = ", vsol)
print("Temperature / K =", temp)

data_up = np.loadtxt('./up/MPF', delimiter='\t')
data_down = np.loadtxt('./down/MPF', delimiter='\t')

data_up = data_up[data_up[:,0].argsort()]
data_down = data_down[data_down[:,0].argsort()]

fig, ax = plt.subplots()

ax.plot(data_up[:,0], data_up[:,1])
ax.plot(data_down[:,0], data_down[:,1])

data = np.zeros(data_up.shape)

for i in range(data_up.shape[0]):
    if(data_up[i,2] > data_down[i,2]):
        data[i,:] = data_down[i,:]
    else:
        data[i,:] = data_up[i,:]


ax.plot(data[:,0], data[:,1])
ax.set_xlabel("mu")
ax.set_ylabel("phi")


VM = np.zeros(data.shape[0])
dlnVM = np.zeros(data.shape[0]-1)

VM = 1.0/data[:,1]*(vsol*6)*1e-27*6.02e23/0.08618 # Molar volume in m3/kg
p = -data[:,2]/(0.25*0.25*0.25) *1e27*8.314/6.02e23/1e6*temp # pressure in MPa

for i in range(VM.shape[0]-1):
    dlnVM[i] = np.log(VM[i])-np.log(VM[i+1])
#    print(i,dlnVM[i])
index = np.argmax(dlnVM)
#print(index,VM[index], VM[index+1], p[index], p[index+1])
p[index+1] = p[index]

print("Pvap / MPa = ", p[index])
print("Vliq / (m^3/kg) = ", VM[index])

fig2,ax2 = plt.subplots()

ax2.plot(VM,p)
ax2.set_ylabel("p / MPa")
ax2.set_xlabel("V / m^3/kg")
ax2.set_xscale('log')
ax2.set_yscale('log')
#plt.show()

plt.show()

Vp =  np.zeros([data.shape[0],2])
Vp[:,0] = VM
Vp[:,1] = p

np.savetxt('Vp', Vp, delimiter='\t')
