import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("RayleighChannelMSEM.csv",delimiter=';')
tACF = df.Theoretical_ACF.tolist()
eACF = df.Empiric_ACF.tolist()
t = df.time.tolist()

k=(10*0.018)/0.0005

plt.plot(t,tACF,'b--',label="Theoretical_ACF")
plt.plot(t,eACF,'g-', label="Empiric_ACF")
plt.grid()
plt.show()
