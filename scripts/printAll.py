import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

# ====================== Rayleigh ======================
ras = pd.read_csv("RayleighSamples.csv",delimiter=';')
rayleigh0 = ras.sig0.to_numpy()
rayleigh1 = ras.sig1.to_numpy()
rayleigh2 = ras.sig2.to_numpy()

ra_x0,ra_y0 = pdf(rayleigh0, B=250)
ra_x1,ra_y1 = pdf(rayleigh1, B=250)
ra_x2,ra_y2 = pdf(rayleigh2, B=250)

# ====================== Rice ======================
ris = pd.read_csv("RiceSamples.csv",delimiter=';')
rice0 = ris.rho0.to_numpy()
rice1 = ris.rho1.to_numpy()
rice2 = ris.rho2.to_numpy()

ri_x0,ri_y0 = pdf(rice0, B=250)
ri_x1,ri_y1 = pdf(rice1, B=250)
ri_x2,ri_y2 = pdf(rice2, B=250)

# ====================== Lognormal ======================
ls = pd.read_csv("LognormalSamples.csv",delimiter=';')
lognormal0 = ls.sig0.to_numpy()
lognormal1 = ls.sig1.to_numpy()
lognormal2 = ls.sig2.to_numpy()

l_x0,l_y0 = pdf(lognormal0)
l_x1,l_y1 = pdf(lognormal1, B=50000)
l_x2,l_y2 = pdf(lognormal2, B=150000)

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

# ====================== Plots ======================

plt.figure()
plt.title("Gaussian Samples")
plt.xlabel("x")
plt.ylabel("PDF")
plt.plot(g_x0,g_y0, label="$\sigma^2 = 1$")
plt.plot(g_x1,g_y1, label="$\sigma^2 = 2$")
plt.plot(g_x2,g_y2, label="$\sigma^2 = 3$")
plt.legend()
plt.grid(linestyle=":")

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
