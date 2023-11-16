import numpy as np
from scipy.integrate import odeint
from scipy import misc, special, integrate
from scipy.interpolate import interp1d
import PlotSettings as ps
#import LatexSettings
from PlotSettings import plt as plt


"""
Constants and Parameters ////////////////////
"""

Mp = 1.22E19 # Planck mass in GeV
me = 0.5e-3 # Electron Mass
mm = 0.105 # Muon Mass
DensityFactor = 1.5E8 # Product of s0/rho_c h^2

# From g-2 mu
gp = 3e-4
mZp = 0.02 # GeV

## Charges

Ql = 1
Qc = 1


# Local directory paths in windows and ubuntu
Thesis_path = 'D:/Nicolás/Documents/Research/Masters-Thesis/Thesis Repo/Thesis-Repo-Public/'
#Thesis_path = '/home/nicolas/Documents/Research/Masters-Thesis/Thesis-Repo-Public/'

"""
Functions ////////////////////
"""

#Annihilation cross-section to leptons
def sigma_to_LL(s, g, M, m):
    Prefac = 4*np.pi*Qc*Qc*Ql
    alpha = g**2/(4*np.pi)
    return Prefac*alpha**2*s*(1 + 2*m**2/s)/((s - M**2)**2*np.sqrt(1 - 4*m**2/s))

#Non relativistic thermal cross-section
def sigvNR(g, M, m):
    return (Qc**4*g**4/np.pi)*(m**2/(4*m**2 - M**2)**2)

#Relativistic d.o.f.
def gs_to_half(T):
    a, b, c = 10.2, 2.349, 0.252
    return a/(1 + np.exp( -b*(T - c)))

#Yield at equilibrium
def Yeq(x):
    return 0.145*x**1.5*np.exp(-x)



x =np.logspace(np.log10(1),np.log10(100), int(1e5))

#Boltzmann equation
def dYdx(Y, x, σv, m):
    return - np.sqrt(np.pi*Mp**2/45) * gs_to_half(m/x)/x**2 * (Y**2-Yeq(x)**2) * m * σv

# Comoving abundance
def Y(x, g, M, m):
    σv = sigvNR(g, M, m)
    return odeint(dYdx, Yeq(x[0]), x, args=(σv, m))


"""
Plots /////////////////////////////
"""

m_array_0 = np.logspace(-3, 1, 1000)

def Omega_h2(g, M, m):
    i_p = 50
    m_array = np.logspace(-3, 3, i_p)
    Y_array = np.zeros(len(m_array))
    for i in range(len(m_array)):
        Y_array[i] = Y(x, g, M, m_array[i])[-1]

    Y_infinity = interp1d(m_array, Y_array, kind='cubic')
    return DensityFactor*Y_infinity(m)*m


g_list = [1e-5, 5e-4, 1e-4, 3e-3]
M_list = [3e-3, 2e-2, 1e-1, 1]
color_list=[ps.MainColor1, ps.BackgroundColor1, ps.MainColor2, ps.Gray1]

fig, ax = plt.subplots(1, 2, figsize=(10,6))



#Left Plot
ax[0].set_title(r'$m_{Z^{\prime}} = 20 ~\mathrm{{MeV}}$')
ax[0].loglog(m_array_0, Omega_h2(g_list[0], M_list[1], m_array_0), label=r'$g^{\prime} = 1\times 10^{-5}$', color=color_list[3])
ax[0].loglog(m_array_0, Omega_h2(g_list[1], M_list[1], m_array_0), label=r'$g^{\prime} = 5\times 10^{-4}$', color=color_list[1])
ax[0].loglog(m_array_0, Omega_h2(g_list[2], M_list[1], m_array_0), label=r'$g^{\prime} = 1\times 10^{-4}$', color=color_list[2])
ax[0].loglog(m_array_0, Omega_h2(g_list[3], M_list[1], m_array_0), label=r'$g^{\prime} = 3\times 10^{-3}$', color=color_list[0])

#Right Plot
ax[1].set_title(r'$g^{\prime} = 5\times 10^{-4}$')
ax[1].loglog(m_array_0, Omega_h2(g_list[1], M_list[0], m_array_0), label=r'$m_{Z^{\prime}} = 3 ~\mathrm{{MeV}}$', color=color_list[0])
ax[1].loglog(m_array_0, Omega_h2(g_list[1], M_list[1], m_array_0), label=r'$m_{Z^{\prime}} = 20 ~\mathrm{{MeV}}$', color=color_list[1])
ax[1].loglog(m_array_0, Omega_h2(g_list[1], M_list[2], m_array_0), label=r'$m_{Z^{\prime}} = 100 ~\mathrm{{MeV}}$', color=color_list[2])
ax[1].loglog(m_array_0, Omega_h2(g_list[1], M_list[3], m_array_0), label=r'$m_{Z^{\prime}} = 1 ~\mathrm{{GeV}}$', color=color_list[3])


for fig_index in range(2):
    ax[fig_index].set_xlabel(r'$m_{\chi} ~[\mathrm{GeV}]$')
    ax[fig_index].set_ylabel(r'$\Omega h^2$')
    ax[fig_index].set_ylim(top=1e5)
    ax[fig_index].axhline(y = 0.12, color='k')
    ax[fig_index].legend(fontsize=14)

fig.tight_layout()
fig.show()
fig.savefig(Thesis_path+'RelicDensity-new.svg')
