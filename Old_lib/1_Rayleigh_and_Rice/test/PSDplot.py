import matplotlib.pyplot as plt
import pandas as pd 

data = pd.read_csv('PSD_test.csv',delimiter=';')
freq = data.frequency.to_list()
JakesPSD = data.JakesPSD.to_list()
GaussianPSD = data.GaussianPSD.to_list()

fig, (ax0,ax1) = plt.subplots(2,1)
ax0.plot(freq,JakesPSD)
ax0.set_title("JakesPSD")
ax0.set_ylabel("PSD(f)")

ax1.plot(freq,GaussianPSD)
ax1.set_title("GaussianPSD")
ax1.set_ylabel("PSD(f)")
ax1.set_xlabel("frequency in Hz")

plt.show()
