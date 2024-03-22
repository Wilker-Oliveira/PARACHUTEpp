import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Lognormal ======================
ls = pd.read_csv("LognormalSamples.csv",delimiter=';')
lognormal0 = ls.sig0.to_numpy()
lognormal1 = ls.sig1.to_numpy()
lognormal2 = ls.sig2.to_numpy()

l_x0,l_y0 = pdf(lognormal0)
l_x1,l_y1 = pdf(lognormal1, B=50000)
l_x2,l_y2 = pdf(lognormal2, B=150000)

plt.figure()
plt.title("Lognormal Samples")
plt.xlabel("x")
plt.ylabel("PDF")
plt.plot(l_x0,l_y0, label="$\sigma^2 = 1$")
plt.plot(l_x1,l_y1, label="$\sigma^2 = 2$")
plt.plot(l_x2,l_y2, label="$\sigma^2 = 3$")
plt.ylim(0,1.2)
plt.xlim(0,6)
plt.legend()
plt.grid(linestyle=":")

plt.show(block = True)
