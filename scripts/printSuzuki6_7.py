import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

def pdf(x, B=3000):
    count, bins_count = np.histogram(x, bins=B, density=True)
    bins_count=np.delete(bins_count,0)
    return bins_count, count

# ====================== Suzuki ======================


# Figure 6.7
s67 = pd.read_csv("SuzukiSamples_figure6_7.csv",delimiter=';')
suzuki6700 = s67.rho0_sigmu0.to_numpy()
suzuki6701 = s67.rho0_sigmu1.to_numpy()
suzuki6702 = s67.rho0_sigmu2.to_numpy()

suzuki6710 = s67.rho1_sigmu0.to_numpy()
suzuki6711 = s67.rho1_sigmu1.to_numpy()
suzuki6712 = s67.rho1_sigmu2.to_numpy()

s670_x0, s670_y0 = pdf(suzuki6700, B=250)
s670_x1, s670_y1 = pdf(suzuki6701, B=700)
s670_x2, s670_y2 = pdf(suzuki6702, B=7000)

s671_x0, s671_y0 = pdf(suzuki6710, B=250)
s671_x1, s671_y1 = pdf(suzuki6711, B=700)
s671_x2, s671_y2 = pdf(suzuki6712, B=7000)

fig, ax1 = plt.subplots(1,1)

plt.title("Suzuki Samples")
plt.xlabel("x")
plt.ylabel("PDF")
p1 = ax1.plot(s670_x0,s670_y0,'-b')
p2 = ax1.plot(s670_x1,s670_y1,'--b')
p3 = ax1.plot(s670_x2,s670_y2,':b')
p4 = ax1.plot(s671_x0,s671_y0,'-m')
p5 = ax1.plot(s671_x1,s671_y1,'--m')
p6 = ax1.plot(s671_x2,s671_y2,':m')
ax1.set_xlim(0,5)

b_patch = Line2D([0],[0],color='b',lw=2, label="$\\rho = 0$")
m_patch = Line2D([0],[0],color='m',lw=2, label="$\\rho = 3$")
sig0_patch = Line2D([0],[0],color='k',linestyle='-', label="$\\sigma_{\mu} = 1/3$")
sig1_patch = Line2D([0],[0],color='k',linestyle='--', label="$\\sigma_{\mu} = 2/3$")
sig2_patch = Line2D([0],[0],color='k',linestyle=':', label="$\\sigma_{\mu} = 1$")
ax1.legend(handles=[b_patch,m_patch,sig0_patch,sig1_patch,sig2_patch])
plt.grid(linestyle=":")

plt.show(block = True)
