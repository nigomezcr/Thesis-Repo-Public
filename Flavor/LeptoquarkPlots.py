import numpy as np
import PlotSettings as ps
from PlotSettings import plt as plt
import LatexSettings as ls

"""
Constants and Patameters /////////////////////////////////
"""
# Constants
GFermi = 1.166e-5 # GeV^-2
alpha = 1/137
V_st = 0.0388
V_bt = 1.013
V_bc = 0.041

# Mass of the LQs
mS3 = 1.2e3
mS1 = 2e3

# Wilson coefficients from fits, along with st. devs

Cmu_fit = 0.67
sigmu_fit = 0.15

Ce_fit = 0.28
sige_fit = 0.2

Cv_fit = 0.11
sigv_fit = 0.02

Cs_fit = 0.06
sigs_fit = 0.05


"""
Yukawa functions /////////////////////////////////
"""

def Yukawa_mu_s(m, y, s):
    Prefactor = np.sqrt(2)*GFermi*alpha*V_st*V_bt/np.pi
    return Prefactor*(Cmu_fit + s)*m**2/y

def Yukawa_e_s(m, y, s):
    Prefactor = np.sqrt(2)*GFermi*alpha*V_st*V_bt/np.pi
    return Prefactor*(Ce_fit + s)*m**2/y

def Yukawa_V_c(m, y, s):
    Prefactor = np.sqrt(2)*GFermi*V_bc
    return Prefactor*(Cv_fit + s)*m**2/y

def Yukawa_S_c(m, y, s):
    Prefactor = np.sqrt(2)*GFermi*V_bc
    return Prefactor*(Cs_fit + s)*m**2/y

"""
Yukawa arrays for the plots  /////////////////////////////////
"""

# For Cmumu plot

yb = np.logspace(-3, 0, 100)
ys = Yukawa_mu_s(mS3, yb, 0)
ys_plus = Yukawa_mu_s(mS3, yb, sigmu_fit)
ys_minus = Yukawa_mu_s(mS3, yb, -sigmu_fit)
ys_plus2 = Yukawa_mu_s(mS3, yb, 2*sigmu_fit)
ys_minus2 = Yukawa_mu_s(mS3, yb, -2*sigmu_fit)


# For Cee plot

yes = Yukawa_e_s(mS3, yb, 0)
yes_plus = Yukawa_e_s(mS3, yb, sige_fit)
yes_minus = Yukawa_e_s(mS3, yb, -sige_fit)
yes_plus2 = Yukawa_e_s(mS3, yb, 2*sige_fit)
yes_minus2 = Yukawa_e_s(mS3, yb, -2*sige_fit)


# For Cv plot

yt = np.linspace(0.2, 1.6, 100)
yc = Yukawa_V_c(mS1, yt, 0)
yc_plus = Yukawa_V_c(mS1, yt, sigv_fit)
yc_minus = Yukawa_V_c(mS1, yt, -sigv_fit)
yc_plus2 = Yukawa_V_c(mS1, yt, 2*sigv_fit)
yc_minus2 = Yukawa_V_c(mS1, yt, -2*sigv_fit)


# For Cs plot

ycp = Yukawa_S_c(mS1, yt, 0)
ycp_plus = Yukawa_S_c(mS1, yt, sigs_fit)
ycp_minus = Yukawa_S_c(mS1, yt, -sigs_fit)
ycp_plus2 = Yukawa_S_c(mS1, yt, 2*sigs_fit)
ycp_minus2 = Yukawa_S_c(mS1, yt, -2*sigs_fit)


"""
Importing data on constrains //////////////////////////
"""
# B -> K nu nu
BtoKnunu_data = np.loadtxt("Data-Sets/BtoKnunu.csv", delimiter=',')
Knu_x = BtoKnunu_data[:,0]
Knu_y = BtoKnunu_data[:,1]

# bbbar -> mu mubar
bbmumu_data = np.loadtxt('Data-Sets/bbmumu.csv', delimiter=',')
bmu_x = bbmumu_data[:,0]
bmu_y = bbmumu_data[:,1]

# B -> tau nu
Btotaunu_data = np.loadtxt("Data-Sets/Btotaunu.csv", delimiter=',')
Tanu_x = Btotaunu_data[:,0]
Tanu_y = Btotaunu_data[:,1]

# ttbar -> tau taubar
tttata_data = np.loadtxt('Data-Sets/ttbartautaubar.csv', delimiter=',')
tta_x = tttata_data[:,0]
tta_y = tttata_data[:,1]

"""
Plotting Functions ///////////////////////////////////
"""
# Function for the plot


def Which_Plot(Set_option):

    Leptoquark, Coefficient = Set_option

    # Preferred regions

    sig2Color = ps.BackgroundColor1
    sig1Color = ps.MainColor1

    sig2Color2 = ps.BackgroundColor2
    sig1Color2 = ps.MainColor2

    # Exclusion region
    if Leptoquark == 'S3':
        DDcolor = ps.Gray1
        ysMin = 0.43
        plt.hlines(ysMin, yb[0], yb[-1], color=  DDcolor)
        plt.fill_between(yb, ysMin, 1, color=DDcolor, alpha=ps.RegionAlpha)
        plt.text(1e-2, 6e-1, '$D - \overline{D}$', color='k')

        KnuColor = ps.Gray2
        plt.plot(Knu_x, Knu_y, color=KnuColor)
        plt.fill_between(Knu_x, Knu_y, 1, color=KnuColor, alpha= ps.RegionAlpha)
        plt.text(3e-1, 2.5e-1, r'$B  \rightarrow K^* \nu \nu $', color='k', rotation=-47 )

        bmuColor = ps.Gray2
        plt.plot(bmu_x, bmu_y, color=bmuColor)
        plt.fill_between(bmu_x, 1e-2, bmu_y, color=bmuColor, alpha= ps.RegionAlpha)
        plt.text(1e-1, 3e-2, r'$b\overline{b}\rightarrow \mu \overline{\mu} $', color='k', rotation= 48)

        pptommColor = ps.Gray1
        ybMin = 0.72
        plt.vlines(ybMin, ys[0], ys[-1], color = pptommColor)
        plt.fill_betweenx(ys, ybMin, 1, color = pptommColor, alpha=ps.RegionAlpha)
        plt.text(7.5e-1, 4e-2, r'$pp \rightarrow \mu \mu $', color='k', rotation='vertical' )


    if Leptoquark == 'S1':
        ppttColor = ps.Gray1
        ycMin = 1
        plt.hlines(ycMin, yt[0], yt[-1], color= ppttColor )
        plt.fill_between(yt, ycMin, 1.6, color=ppttColor, alpha=ps.RegionAlpha)
        plt.text(0.3, 1.3, r'$pp \rightarrow \tau \tau$',  color='k')

        ZtautauColor = ps.Gray2
        ybMin = 1.2
        plt.vlines(ybMin, yt[0], yt[-1], color = ZtautauColor)
        plt.fill_betweenx(yt, ybMin, 1.6, color = ZtautauColor, alpha=ps.RegionAlpha)
        plt.text(1.22, 0.5, r'$Z \rightarrow \tau \tau $',  color= 'k', rotation='vertical' )


        TaunuColor = ps.Gray1
        plt.plot(Tanu_x, Tanu_y, color=TaunuColor)
        plt.fill_between(Tanu_x, Tanu_y, 1.6, color=TaunuColor, alpha= ps.RegionAlpha)
        plt.text(0.7, 1.2, r'$B  \rightarrow K^* \nu \nu $',  color='k', rotation=-48 )


        tttauColor = ps.Gray2
        plt.plot(tta_x, tta_y, color=tttauColor)
        plt.fill_between(tta_x, 0.2, tta_y, color=tttauColor, alpha= ps.RegionAlpha)
        plt.text(0.7, 0.85, r'$ tt \rightarrow \tau \tau $',  color='k', rotation= 45)


    if Coefficient == 'Cee':
        plt.title('$C_9^{ee} = - C_{10}^{ee}$')
        plt.fill_between(yb, yes_minus2, yes_plus2, color=sig2Color2)
        plt.fill_between(yb, yes_minus, yes_plus, color=sig1Color2)

        plt.plot(yb, yes, color=sig1Color2)
        plt.text(1.5e-3, 2.5e-2, r'$m_S = {} $ TeV'.format(mS3/1e3), color='black')
        patch1 = ps.mpatches.Patch(color=sig1Color2, label='$R_{K^*} (1 \sigma)$')
        patch2 = ps.mpatches.Patch(color=sig2Color2, label='$R_{K^*}  (2 \sigma)$')

        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-2, 1)
        plt.xlim(1e-3, 1)
        plt.xlabel('$-y_b$')
        plt.ylabel('$y_s$')


    if Coefficient == 'Cmm':
        plt.title('$C_9^{\mu \mu} = - C_{10}^{\mu \mu}$')
        plt.fill_between(yb, ys_minus2, ys_plus2, color=sig2Color2)
        plt.fill_between(yb, ys_minus, ys_plus, color=sig1Color2)

        plt.plot(yb, ys, color=sig1Color2)
        plt.text(1.5e-3, 2.5e-2, r'$m_S = {} $ TeV'.format(mS3/1e3), color='black')
        patch1 = ps.mpatches.Patch(color=sig1Color2, label='$R_{K^*} (1 \sigma)$')
        patch2 = ps.mpatches.Patch(color=sig2Color2, label='$R_{K^*}  (2 \sigma)$')

        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-2, 1)
        plt.xlim(1e-3, 1)
        plt.xlabel('$-y_b$')
        plt.ylabel('$y_s$')

    if Coefficient == 'Cs':
        plt.title('$C_s = - 4 C_{T}$')
        plt.fill_between(yt, ycp_minus2, ycp_plus2, color=sig2Color)
        plt.fill_between(yt, ycp_minus, ycp_plus, color=sig1Color)

        plt.text(0.25, 0.5, r'$m_S = {} $ TeV'.format(mS1/1e3), color='black')
        patch1 = ps.mpatches.Patch(color=sig1Color, label='$R_{D^{(*)}} (1 \sigma)$')
        patch2 = ps.mpatches.Patch(color=sig2Color, label='$R_{D^{(*)}}  (2 \sigma)$')

        plt.xscale('linear')
        plt.yscale('linear')
        plt.ylim(0.2, 1.6)
        plt.xlim(0.2, 1.4)
        plt.xlabel(r'$y_{\tau b}$')
        plt.ylabel(r'$y_c$')


    if Coefficient == 'Cv':
        plt.title('$C_{V_L}$')
        plt.fill_between(yt, yc_minus2, yc_plus2, color=sig2Color)
        plt.fill_between(yt, yc_minus, yc_plus, color=sig1Color)

        plt.text(0.25, 0.5, r'$m_S = {} $ TeV'.format(mS1/1e3), color='black')
        patch1 = ps.mpatches.Patch(color=sig1Color, label='$R_{D^(*)} (1 \sigma)$')
        patch2 = ps.mpatches.Patch(color=sig2Color, label='$R_{D^(*)}  (2 \sigma)$')

        plt.xscale('linear')
        plt.yscale('linear')
        plt.ylim(0.2, 1.6)
        plt.xlim(0.2, 1.4)
        plt.xlabel(r'$y_{\tau b}$')
        plt.ylabel('$y\'_c$')




    plt.legend(handles=[patch1, patch2], loc='lower left')

# Define the plots to make
Option1 = ('S3', 'Cmm')
Option2 = ('S3', 'Cee')
Option3 = ('S1', 'Cv')
Option4 = ('S1', 'Cs')

# Save the plots
def Save_Plot(Set_option, Format):

    if Set_option == Option1:
        Which_Plot(Set_option)
        plt.savefig('FlavourPlots/LQ-S3-Cmm.'+Format)

    if Set_option == Option2:
        Which_Plot(Set_option)
        plt.savefig('FlavourPlots/LQ-S3-Cee.'+Format)

    if Set_option == Option3:
        Which_Plot(Set_option)
        plt.savefig('FlavourPlots/LQ-S1-Cv.'+Format)

    if Set_option == Option4:
        Which_Plot(Set_option)
        plt.savefig('FlavourPlots/LQ-S1-Cs.'+Format)


# Output final functions
pdf = 'pdf' # Format for plot
svg = 'svg'

Save_Plot(Option1, svg)
#Save_Plot(Option2, svg)
#Save_Plot(Option3, svg)
#Save_Plot(Option4, svg)
