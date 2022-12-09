import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift
from mpl_toolkits.mplot3d import Axes3D

#usefull functions
def lin2db(x):
    return 10.0*np.log10(x)

def db2lin(x):
    return 10.0**(x/10.0)

def lin2dbm(x):
    return 10.0*np.log10(x)+30.0

def dbm2lin(x):
    return 10.0**(x/10.0 - 3.0)

def AutoCorrelation(h,k):
    k=int(k)
    ri=np.zeros(k,dtype = 'complex_')
    for i in range(k):
        for ki in range(k-i):
            ri[i]+=h[ki+i]*h[ki].conj().T

    return ri

df=pd.read_csv("RayleighChannelMC.csv",delimiter=';')

u1 = np.array(df.u1.tolist())
u2 = 1j*np.array(df.u2.tolist())
t = np.array(df.time.tolist())

coherence_time=float(input("coherence time: "))

k=(10*coherence_time)/0.0005


rij=np.zeros(int(k/2),dtype = 'complex_')

h=u1+u2

g=(np.abs(h)**2)

print(np.mean(g))
print(np.mean(np.imag(h)))
print(np.mean(np.real(h)))

plt.figure()
plt.plot(t[0:int(k)],lin2db(g[0:int(k)]))

plt.figure()
rij=AutoCorrelation(h,k/2)
plt.plot(np.real(rij)/len(rij))

plt.show()
