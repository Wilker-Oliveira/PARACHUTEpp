import matplotlib.pyplot as plt
import pandas as pd 

data = pd.read_csv('ACF_test.csv',delimiter=';')
tau = data.tau.to_list()
JakesACF = data.JakesACF.to_list()
GaussianACF = data.GaussianACF.to_list()

fig, (ax0,ax1) = plt.subplots(2,1)
ax0.plot(tau,JakesACF)
ax0.set_title("JakesACF")
ax0.set_ylabel("ACF(tau)")

ax1.plot(tau,GaussianACF)
ax1.set_title("GaussianACF")
ax1.set_ylabel("ACF(tau)")
ax1.set_xlabel("tau in s")

plt.show()
