import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Rice ======================
ris = pd.read_csv("RiceSamples.csv",delimiter=';')
rice0 = ris.rho0.to_numpy()
rice1 = ris.rho1.to_numpy()
rice2 = ris.rho2.to_numpy()

ri_x0,ri_y0 = pdf(rice0, B=250)
ri_x1,ri_y1 = pdf(rice1, B=250)
ri_x2,ri_y2 = pdf(rice2, B=250)

plt.figure()
plt.title("Rice Samples")
plt.xlabel("x")
plt.ylabel("PDF")
plt.plot(ri_x0,ri_y0, label="$\\rho=0$")
plt.plot(ri_x1,ri_y1, label="$\\rho=2$")
plt.plot(ri_x2,ri_y2, label="$\\rho=3$")
plt.ylim(0,0.7)
plt.xlim(0,6)
plt.legend()
plt.grid(linestyle=":")

plt.show(block = True)
