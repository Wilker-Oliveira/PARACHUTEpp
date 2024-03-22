import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Rayleigh ======================
ras = pd.read_csv("RayleighSamples.csv",delimiter=';')
rayleigh0 = ras.sig0.to_numpy()
rayleigh1 = ras.sig1.to_numpy()
rayleigh2 = ras.sig2.to_numpy()

ra_x0,ra_y0 = pdf(rayleigh0, B=250)
ra_x1,ra_y1 = pdf(rayleigh1, B=250)
ra_x2,ra_y2 = pdf(rayleigh2, B=250)

plt.figure()
plt.title("Rayleigh Samples")
plt.xlabel("x")
plt.ylabel("PDF")
plt.plot(ra_x0,ra_y0, label="$\sigma^2 = 1$")
plt.plot(ra_x1,ra_y1, label="$\sigma^2 = 2$")
plt.plot(ra_x2,ra_y2, label="$\sigma^2 = 3$")
plt.ylim(0,0.7)
plt.xlim(0,6)
plt.legend()
plt.grid(linestyle=":")

plt.show(block = True)
