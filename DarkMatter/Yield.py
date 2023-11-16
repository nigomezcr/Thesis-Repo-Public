import numpy as np
from scipy.integrate import odeint
import PlotSettings as ps
#import LatexSettings
from PlotSettings import plt as plt
from Constants import Mp
from CrossSections import sigvNR
"""
Constants and Parameters ////////////////////
"""
# From g-2 mu
gp = 3e-4
mZp = 0.02 # GeV


# Local directory paths in windows and ubuntu
Thesis_path = 'D:/Nicolás/Documents/Research/Masters-Thesis/Thesis Repo/Thesis-Repo/'
#Thesis_path = '/home/nicolas/Documents/Research/Masters-Thesis/Thesis-Repo/'

"""
//////////////////// Yield Functions ////////////////////
"""
#Relativistic d.o.f.
def gs_to_half(T):
    a, b, c = 10.2, 2.349, 0.252
    return a/(1 + np.exp( -b*(T - c)))

#Yield at equilibrium
def Yeq(x):
    return 0.145*x**1.5*np.exp(-x)



x =np.logspace(np.log10(1),np.log10(100), int(1e3))

#Boltzmann equation
def dYdx(Y, x, σv, m):
    return - np.sqrt(np.pi*Mp**2/45) * gs_to_half(m/x)/x**2 * (Y**2-Yeq(x)**2) * m * σv

# Comoving abundance
def Y(x, g, M, m):
    σv = sigvNR(g, M, m)
    return odeint(dYdx, Yeq(x[0]), x, args=(σv, m))


"""
 //////////////////////////////// Plots Functions //////////////////////////////
"""
plt.loglog(x,Yeq(x), color='k')

M = 5e-2

m = 1
gp = 3e-3
Y1 = Y(x, gp, M, m)
plt.loglog(x, Y1, label=r'$m_{{\chi}} = {} ~\mathrm{{GeV}},g^{{\prime}} = {} \times 10^{{-3}} $'.format(m, gp*1e3), color=ps.MainColor1)
plt.text(22, 1.4e-2, r'$\sigma v= 1.6 \times 10^{-12} ~\mathrm{GeV^{-2}}$', color=ps.Gray1 , size=10)

m = 100
gp = 3e-3
Y1 = Y(x, gp, M, m)
plt.loglog(x, Y1, label=r'$m_{{\chi}} = {} ~\mathrm{{GeV}},g^{{\prime}} = {}  \times 10^{{-3}} $'.format(m, gp*1e3), color=ps.BackgroundColor2)

gp = 5e-4
Y1 = Y(x, gp, M, m)
plt.loglog(x, Y1, label=r'$m_{{\chi}} = {} ~\mathrm{{GeV}},g^{{\prime}} = {}  \times 10^{{-3}}$'.format(m, gp*1e3), color=ps.BackgroundColor1)

m = 1
gp = 5e-4
Y1 = Y(x, gp, M, m)
plt.loglog(x, Y1, label=r'$m_{{\chi}} = {} ~\mathrm{{GeV}},g^{{\prime}} = {}  \times 10^{{-3}}$'.format(m, gp*1e3), color=ps.MainColor2)
plt.text(22, 2.5e-6, r'$\sigma v= 1.2 \times 10^{-15} ~\mathrm{GeV^{-2}}$', color=ps.Gray1, size=10)

plt.text(1.2, 1e-9, r'$m_{Z^{\prime}} = 50 ~\mathrm{MeV}$')

plt.ylim(1E-14,1)
plt.xlim(1,100)
plt.xlabel('$x = m_{\chi}/T$')
plt.ylabel('$Y$')
plt.legend(loc=3, fontsize=12)
#plt.savefig(Thesis_path+'Yield.svg')


print('Done')
