import numpy as np
import matplotlib.pyplot as plt
import PlotSettings as ps
plt.rcdefaults()

# Anomalous magnetic moment results from theory and experiment
aSM = 18.04
sig_aSM = 0.51

aBNL = 20.8
sig_aBNL = 0.63

aFNAL21 = 20.40
sig_aFNAL21 = 0.54

aFNAL23 = 20.55
sig_aFNAL23 = 0.24

aEXP = 20.59
sig_aEXP = 0.22

# Colors for the plots
SM_color = ps.MainColor1
BNL_color = ps.Gray1
FNAL_color = ps.MainColor2
EXP_color = ps.BackgroundColor2
alph = ps.RegionAlpha


#Settings for the plots
font = {'family': 'dosis',
        'color':  'darkred',
        'weight': 'normal',
        'size': 12,
        }


plt.xlim(xmin=17, xmax=22)
plt.ylim(ymin=0.8, ymax=2)
plt.xlabel(r'$a_{\mu} \times 10^9 - 1165900$')
plt.tick_params(left = False, right = False , labelleft = False)
plt.tight_layout()

#Functions for the plots

#SM
def plot_SM():
    plt.plot(aSM, 1, 'o', color=SM_color)
    plt.hlines(1, xmin=aSM - sig_aSM, xmax=aSM + sig_aSM, colors=SM_color)
    plt.fill_betweenx(y=[0.8, 2], x1=aSM - sig_aSM, x2=aSM + sig_aSM, color=SM_color, alpha=alph)
    plt.text(aSM - sig_aSM, 0.9, 'Standard Model', color=SM_color, fontdict=font)

#BNL
def plot_BNL():
    plt.plot(aBNL, 1.8, 'v', color=BNL_color)
    plt.hlines(1.8, xmin=aBNL - sig_aBNL, xmax=aBNL + sig_aBNL, colors=BNL_color)
    plt.text(aBNL - 0.2, 1.85, 'BNL', color=BNL_color, fontdict=font)
    plt.savefig('FlavourPlots/MuonMoment-Measurements1.svg')

#FNAl '21
def plot_FNAL21():
    plt.plot(aFNAL21, 1.6, 'o', color=FNAL_color)
    plt.hlines(1.6, xmin=aFNAL21 - sig_aFNAL21, xmax=aFNAL21 + sig_aFNAL21, colors=FNAL_color)
    plt.text(aFNAL21 - sig_aFNAL21, 1.65, 'FNAL Run 1', color=FNAL_color, fontdict=font)
    plt.savefig('FlavourPlots/MuonMoment-Measurements2.svg')


#FNAL '23
def plot_FNAL23():
    plt.plot(aFNAL23, 1.4, 'o', color=FNAL_color)
    plt.hlines(1.4, xmin=aFNAL23 - sig_aFNAL23, xmax=aFNAL23 + sig_aFNAL23, colors=FNAL_color)
    plt.text(aFNAL23 - sig_aFNAL23, 1.45, 'FNAL Run 2/3', color=FNAL_color, fontdict=font)
    plt.savefig('FlavourPlots/MuonMoment-Measurements3.svg')

# EXP
def plot_EXP_average():
    plt.plot(aEXP, 1.2, 'D', color=EXP_color)
    plt.hlines(1.2, xmin=aEXP - sig_aEXP, xmax=aEXP + sig_aEXP, colors=EXP_color)
    plt.fill_betweenx(y=[0.8, 2], x1=aEXP - sig_aEXP, x2=aEXP + sig_aEXP, color=EXP_color, alpha=alph)
    plt.text(aEXP - sig_aEXP - 0.2, 1.25, 'Exp Average', color=EXP_color, fontdict=font)
    plt.savefig('FlavourPlots/MuonMoment-Measurements.svg')


def Plot_All():
    plot_SM()
    plot_BNL()
    plot_FNAL21()
    plot_FNAL23()
    plot_EXP_average()

Plot_All()
