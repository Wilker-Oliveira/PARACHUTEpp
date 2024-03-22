import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Suzuki ======================

# Figure 2.4
s24 = pd.read_csv("SuzukiSamples_figure2_4.csv",delimiter=';')
suzuki2400 = s24.sig0_sigmu0.to_numpy()
suzuki2401 = s24.sig0_sigmu1.to_numpy()
suzuki2402 = s24.sig0_sigmu2.to_numpy()

suzuki2410 = s24.sig1_sigmu0.to_numpy()
suzuki2411 = s24.sig1_sigmu1.to_numpy()
suzuki2412 = s24.sig1_sigmu2.to_numpy()

s240_x0, s240_y0 = pdf(suzuki2400, B=250)
s240_x1, s240_y1 = pdf(suzuki2401, B=250)
s240_x2, s240_y2 = pdf(suzuki2402, B=1200)

s241_x0, s241_y0 = pdf(suzuki2410, B=250)
s241_x1, s241_y1 = pdf(suzuki2411, B=250)
s241_x2, s241_y2 = pdf(suzuki2412, B=1200)

fig, ax1 = plt.subplots(1,1)

plt.title("Suzuki Samples")
plt.xlabel("x")
plt.ylabel("PDF")
p1 = ax1.plot(s240_x0,s240_y0,'-b')
p2 = ax1.plot(s240_x1,s240_y1,'--b')
p3 = ax1.plot(s240_x2,s240_y2,':b')
p4 = ax1.plot(s241_x0,s241_y0,'-m')
p5 = ax1.plot(s241_x1,s241_y1,'--m')
p6 = ax1.plot(s241_x2,s241_y2,':m')
ax1.set_ylim(0,0.7)
ax1.set_xlim(0,10)

b_patch = Line2D([0],[0],color='b',lw=2, label="$\\sigma^2_0 = 1$")
m_patch = Line2D([0],[0],color='m',lw=2, label="$\\sigma^2_0 = 3$")
sig0_patch = Line2D([0],[0],color='k',linestyle='-', label="$\\sigma_{\mu}^2 = 0$")
sig1_patch = Line2D([0],[0],color='k',linestyle='--', label="$\\sigma_{\mu}^2 = 1/3$")
sig2_patch = Line2D([0],[0],color='k',linestyle=':', label="$\\sigma_{\mu}^2 = 1/2$")
ax1.legend(handles=[b_patch,m_patch,sig0_patch,sig1_patch,sig2_patch])
plt.grid(linestyle=":")

plt.show(block = True)
