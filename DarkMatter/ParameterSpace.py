import numpy as np
from scipy.integrate import quad, odeint
from scipy import misc, special, optimize
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd


""" 
Global Constants and parameters /////////////////////////////
"""
# Windows
#Thesis_path = 'D:/Nicolás/Documents/Research/Masters-Thesis/Thesis Repo/Thesis-Repo/'
# Linux
Thesis_path = '/home/nicolas/Documents/Research/Masters-Thesis/Thesis-Repo/'


Mp = 1.22E19 # Planck mass in GeV
me = 0.5e-3 # Electron Mass
mm = 0.105 # Muon Mass
DensityFactor = 1.5e8

#Grid
g_p = 20
M_l = np.logspace(-3, 1, g_p)
g_l = np.logspace(-4, -1, g_p)
g_g, M_g = np.meshgrid(g_l, M_l)


""" 
Functions  //////////////////////////
"""

#### g - 2

#Functions
def integral(m):
    return quad(lambda x: mm**2*x**2*(1-x)/(x**2*mm**2 + (1-x)*m**2), 0, 1)[0]

def Delta_a(g, M, charge):    
    return g**2*charge/(4*np.pi**2) * integral(M)

#Delta a grid
def Delta_grid(charges):
    Da_grid = np.zeros((g_p, g_p))
    for i in range(g_p):
        for j in range(g_p):
            Da_grid[i, j] = Delta_a(g_l[j], M_l[i], charges[0])
    return Da_grid

#### Relic density

def sigvNR(g, M, m, charges):
    return charges[0]*charges[1]*(g**4/np.pi)*(m**2/(4*m**2 - M**2)**2)

def gs_to_half(T):
    a, b, c = 10.2, 2.349, 0.252
    return a/(1 + np.exp( -b*(T - c)))

def Yeq(x):
    return 0.145*x**1.5*np.exp(-x)

x =np.logspace(np.log10(1),np.log10(100), int(1e6))

def dYdx(Y, x, σv, m):
    return - np.sqrt(np.pi*Mp**2/45) * gs_to_half(m/x)/x**2 * (Y**2-Yeq(x)**2) * m * σv

def Y(x, g, M, m, charges):
    σv = sigvNR(g, M, m, charges)
    return odeint(dYdx, Yeq(x[0]), x, args=(σv, m))


def Log_Omega_h2(m, charges):
    Oh2_g = np.zeros((g_p, g_p))
    for i in range(g_p):
        for j in range(g_p):
            Oh2_g[i,j] = Y(x, g_l[j], M_l[i], m, charges)[-1,0]
    return np.log10(Oh2_g*DensityFactor*m)


#for g in g_l:
#    ylist =  Y(x, g, 5e-2, 2e-2, (1,1))[-1, 0]

j = 11
g = g_l[j]
print(g, '\t', Y(x, g, 5e-2, 2e-2, (1,1))[-1,0])

"""
Importing data on constrains //////////////////////////
"""

# BaBar
# BaBar_data = np.loadtxt("Data-Sets/Babar_exclusion.csv", delimiter=',')  # Linux direction
BaBar_data = np.loadtxt("../Flavour/Data-Sets/Babar_exclusion.csv", delimiter=',') # Windows direction
Babar_x = BaBar_data[:,0]
Babar_y = BaBar_data[:,1]

# CCFR
#CCFR_data = np.loadtxt(Thesis_path+"Flavour/Data-Sets/Nu-trident_exclusion.csv", delimiter=',')  # Windows direction
CCFR_data = np.loadtxt("../Flavour/Data-Sets/Nu-trident_exclusion.csv", delimiter=',')  # Linux direction
CCFR_x = CCFR_data[:,0]
CCFR_y = CCFR_data[:,1]


"""
Plotting Functions ///////////////////////////////////
"""

#Colors:
RegionAlpha = 0.5
MainColor1 = (0.580,0.706,0.231)
MainColor2 = (0.650,0.110,0.192)
BackgroundColor1 = (0.275,0.419,0.247)
BackgroundColor2 = (0.467,0.137,0.184)
Gray1 = (0.337,0.352,0.360)
Gray2 = (0.694,0.698,0.690)



"""
#Plots parameters
#Latex rendering
#plt.rc('font', **{'family': 'serif', 'serif': ['Modern Computer']})
#plt.rc('text', usetex=True)
#plt.rc('font', weight='bold')

# Global settings
sns.set_style('ticks') # darkgrid, white grid, dark, white and ticks
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.85 )   # width of the frame
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('font', size=14)          # controls default text sizes

# Parameters for the plots

m_chi = [1e-3, 5e-1, 2e-2, 3]
Qs = [(1, 1), (1, 1), (0.6, 0.8), (0.6, 0.8)] 


#### Parameter space plot

fig, ax = plt.subplots(2,2, figsize=(10,8))
ax = ax.ravel()

delta_a_minus = (259 - 2*51)*1e-11
delta_a_plus = (259 + 2*51)*1e-11

levels_a = [delta_a_minus, delta_a_plus]
levels_d = [np.log10(0.1), np.log10(0.16)]


patch_a = mpatches.Patch(color=MainColor1, label=r'$\Delta a $')
patch_d = mpatches.Patch(color=MainColor2, label=r'$ \Omega h^2 $')


for fig_index in range(4):
    #ax[fig_index].text(2, 2e-4, r'$m_{{\chi}} = {}, Q_{{ \mu}} = {}, \n Q_{{\chi}} = {}$'.format(m_chi[fig_index],Qs[fig_index][0], Qs[fig_index][1]), fontsize=10)
    ax[fig_index].set_xlabel('$m_{Z\'} (\mathrm{GeV} )$' )
    ax[fig_index].set_ylabel('$g\'$')
    ax[fig_index].set_xscale('log')
    ax[fig_index].set_yscale('log')
    ax[fig_index].set_ylim(1e-4, 1e-2)
    ax[fig_index].set_xlim(1e-3, 1e1)

    # Plot exclusion regions
    ax[fig_index].plot(CCFR_x, CCFR_y, color=Gray1, linestyle='--')
    ax[fig_index].vlines(x=5.3e-3, ymin=1e-5, ymax=1e-2, colors=Gray1, linestyle='-')
    ax[fig_index].fill_betweenx(g_l, 5.3e-3, color=Gray2, alpha=RegionAlpha) 
    ax[fig_index].fill_between(CCFR_x, CCFR_y, g_l[-1], color=Gray2, alpha=RegionAlpha) 
    
    #g - 2
    ax[fig_index].contourf(M_g, g_g, Delta_grid(Qs[fig_index]), levels_a, colors= (MainColor1,))

    #Relic Density
    ax[fig_index].contourf(M_g, g_g, Log_Omega_h2(m_chi[fig_index], Qs[fig_index]), levels_d, colors=(MainColor2,))

    ax[fig_index].legend(handles=[patch_a, patch_d], loc='lower right')

#Plot constrains


fig.tight_layout()
fig.show()

"""
