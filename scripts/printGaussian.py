import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Gaussian ======================
gs = pd.read_csv("GaussianSamples.csv",delimiter=';')
gaussian0 = gs.sig0.to_numpy()
gaussian1 = gs.sig1.to_numpy()
gaussian2 = gs.sig2.to_numpy()

g_x0,g_y0 = pdf(gaussian0, B=250)
g_x1,g_y1 = pdf(gaussian1, B=250)
g_x2,g_y2 = pdf(gaussian2, B=250)

plt.figure()
plt.title("Gaussian Samples")
plt.xlabel("x")
plt.ylabel("PDF")
plt.plot(g_x0,g_y0, label="$\sigma^2 = 1$")
plt.plot(g_x1,g_y1, label="$\sigma^2 = 2$")
plt.plot(g_x2,g_y2, label="$\sigma^2 = 3$")
plt.legend()
plt.grid(linestyle=":")

plt.show(block = True)
