import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("RayleighChannelJakes.csv",delimiter=';')
tACF = np.array(df.Theoretical_ACF.tolist())
eACF = np.array(df.Empirical_ACF.tolist())
t = df.time.tolist()

k=(5*0.018)/0.0005

plt.plot(t,tACF/tACF[0],'b--',label="Theoretical_ACF")
plt.plot(t,eACF/eACF[0],'g-', label="Empiric_ACF")
plt.grid()
plt.show()
