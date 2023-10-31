import numpy as np
from scipy.integrate import quad
import PlotSettings as ps
from PlotSettings import plt as plt
import LatexSettings as ls

"""
Constants and Patameters /////////////////////////////////
"""
# Constants
mm = .105 # muon mass in GeV
me = 0.000511 # electron mass in GeV
Delta_amu = 248 # discrepancy in a_mu 10^{-11} units
sigma_amu = 49
Delta_ae = 0.1053 # discrepancy in a_e in 10^{-11} units
sigma_ae = 0.0775

#Parameters for the plots
rho = np.array([0, 0, -1/6, -1/2, -1/20])
nu = np.array([0, 1, 1/3, 1/2, 1/5])

Qe = rho
Qmu = nu - rho


"""
////////// Coupling constant function /////////////////////////////////
"""

def integral(M, m):
    return quad(lambda x: m**2*x**2*(1-x)/(x**2*m**2 + (1-x)*M**2), 0, 1)[0]

def gp(charge, Delta_a, M, m):
    gp2 = 4*np.pi**2*(Delta_a)/( np.abs(charge)*integral(M, m) )
    return np.sqrt(gp2)


"""
Import data on constrains //////////////////////////
"""

# BaBar
BaBar_data = np.loadtxt("Data-Sets/Babar_exclusion.csv", delimiter=',')
Babar_x = BaBar_data[:,0]
Babar_y = BaBar_data[:,1]

# CCFR
CCFR_data = np.loadtxt("Data-Sets/Nu-trident_exclusion.csv", delimiter=',')
CCFR_x = CCFR_data[:,0]
CCFR_y = CCFR_data[:,1]

"""
///////////  Arrays For the Plots  //////////////////////////
"""

N=1000
mZp = np.linspace(1, 7000, N)*10**(-3) # Mass in GeV

#Fill g-2 region for Lmu - Ltau, left plot
g_sigma1p = [gp(Qmu[1], Delta_amu + 1*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
g_sigma1m = [gp(Qmu[1], Delta_amu - 1*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
g_sigma2p = [gp(Qmu[1], Delta_amu + 2*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
g_sigma2m = [gp(Qmu[1], Delta_amu - 2*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]

#Fill g-2 muon array, right plot
g_muon1 = [gp(Qmu[1], Delta_amu, M, mm)*10**(-5.5) for M in mZp]
g_muon2 = [gp(Qmu[2], Delta_amu, M, mm)*10**(-5.5) for M in mZp]
g_muon3 = [gp(Qmu[3], Delta_amu, M, mm)*10**(-5.5) for M in mZp]
g_muon4 = [gp(Qmu[4], Delta_amu, M, mm)*10**(-5.5) for M in mZp]

#Fill g-2 electron array, right plot
#g_electron1 = [gp(Qe[1], Delta_ae, M, me)*10**(-5.5) for M in mZp] No need beacuse Qe=0
g_electron2 = [gp(Qe[2], Delta_ae, M, me)*10**(-5.5) for M in mZp]
g_electron3 = [gp(Qe[3], Delta_ae, M, me)*10**(-5.5) for M in mZp]
g_electron4 = [gp(Qe[4], Delta_ae, M, me)*10**(-5.5) for M in mZp]


"""
Plotting Functions ///////////////////////////////////
"""

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Left Plot
#   CCFR exclusion
CCFRColor = ps.Gray2
axes[0].plot(CCFR_x,CCFR_y,'-',color=CCFRColor, alpha=1)
axes[0].fill_between(CCFR_x, CCFR_y, 1e-1, color=CCFRColor,alpha=ps.RegionAlpha)
axes[0].text(1e-2,3e-2, 'CCFR', fontsize=12, color='k')


# Plot g-2 prefered region
gColor1 = ps.BackgroundColor1
gColor2 = ps.MainColor1
g_patch = ps.mpatches.Patch(color=gColor1, label='$1\sigma$')
g_patch2 = ps.mpatches.Patch(color=gColor2, label='$2\sigma$')

axes[0].set_title(r'$L_{\mu} - L_{\tau}$ model')
axes[0].fill_between(mZp, g_sigma2m, g_sigma2p, color=gColor2)
axes[0].fill_between(mZp, g_sigma1m, g_sigma1p, color=gColor1)
axes[0].legend(handles=[g_patch, g_patch2], loc='lower right')


### Right Plot
#  (g-2)e exclusion

#axes[1].plot(mZp, g_electron1, color='black', linestyle='dotted')
axes[1].plot(mZp, g_electron2, color='black', linestyle='dashed')
axes[1].plot(mZp, g_electron3, color='black', linestyle='dashdot')
axes[1].plot(mZp, g_electron4, color='black', linestyle='solid')
axes[1].fill_between(mZp, g_electron4, 1e-1,color=CCFRColor, alpha=ps.RegionAlpha)
axes[1].text(1e-2 , 3e-2, '$(g-2)_e$', fontsize= 12,color='k')


# Plot g-2 prefered region
axes[1].set_title(r'$ \Lambda$ model with $\nu=1$, $\rho$ varied')
axes[1].plot(mZp, g_muon1, color=gColor1, linestyle='dotted' , label=r'$\rho = 0$')
axes[1].plot(mZp, g_muon2, color=gColor1, linestyle='dashed' , label=r'$\rho = -1/6$')
axes[1].plot(mZp, g_muon3, color=gColor1, linestyle='dashdot', label=r'$\rho = -1/2$')
axes[1].plot(mZp, g_muon4, color=gColor1, linestyle='solid', label=r'$\rho = -1/20$')
axes[1].legend(loc='lower right')


#### Both Plots

for i in (0, 1):
    #Settings
    axes[i].set_ylim(1e-4, 1e-1)
    axes[i].set_xlim(1e-3, mZp[-1])
    axes[i].set_xlabel('$m_{Z\'} ~[\mathrm{GeV}]$')
    axes[i].set_ylabel(' $g\' $ ')
    axes[i].set_yscale('log')
    axes[i].set_xscale('log')

    # Plot joint exclusion regions

    #  BaBar
    BaBarColor = ps.BackgroundColor2
    axes[i].fill_between(Babar_x, Babar_y, 1e-1, color=BaBarColor, alpha=ps.RegionAlpha)
    axes[i].loglog(Babar_x,Babar_y, color=BaBarColor, alpha=0.7)
    axes[i].text(1, 3e-2, 'BaBar', fontsize=12, color=BaBarColor)

    #  Neff
    NeffColor = ps.Gray1
    axes[i].fill_between(([1e-3,5.3e-3]),1e-4,1e-1, color=NeffColor, alpha=ps.RegionAlpha)
    axes[i].text(1.4e-3 , 3e-2, r'$\Delta N_{eff}$', fontsize= 12,color=NeffColor)

fig.tight_layout()

plt.savefig('FlavourPlots/g-2.pdf')
plt.savefig('FlavourPlots/g-2.svg')
