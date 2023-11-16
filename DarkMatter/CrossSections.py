import  numpy as np
from scipy.integrate import quad
from Constants import c, fc, GeVtocm2

"""
///////////// Annihilation Cross Sections  ///////////////
"""

#Annihilation cross-section to leptons
def sigma_to_LL(s, g, M, m):
    Qc = 1
    Ql = 1
    Prefac = 4*np.pi*Qc*Qc*Ql
    alpha = g**2/(4*np.pi)
    return Prefac*alpha**2*s*(1 + 2*m**2/s)/((s - M**2)**2*np.sqrt(1 - 4*m**2/s))

#Non relativistic thermal cross-section
def sigvNR(g, M, m):
    Qc = 1
    return (Qc**4*g**4/np.pi)*(m**2/(4*m**2 - M**2)**2)


"""
///////////// Scattering Cross Sections Computed  ///////////////
"""


#### Total cross section, analitical #####
#  Particle-Particle: chi chi -> chi chi
def sigma_r(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = (2*g**4/(d*np.pi*s*m))
    Term1 = ((2*s + 3*M**2)*s*v**2 + 2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2 + s*v**2))
    LogTerm1 = -(s*v**2*(3*M**2 + 4*m**2) + 2*(M**2 + 2*m**2)**2 - 4*m**4)/(s*v**2*(2*M**2 + s*v**2))
    return  fc*Prefactor*(Term1 + LogTerm1*np.log(1 + s*v**2/M**2) )


# Particle-Antiparticle: chi chibar -> chi chibar
def sigma_a(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = g**4/(d*np.pi*s*m)
    sigma_t = Prefactor*( (s*v**2*(2*s + 3*M**2)+2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2+s*v**2)) - (M**2 + s)/(s*v**2)*np.log(1 + s*v**2/M**2) )
    sigma_s = Prefactor/3*( (12*m**4 + 6*(v**2 + 1)*s*m**2 + s**2*v**4 ) / (s - M**2)**2 )
    sigma_st = -Prefactor/2*( (16*m**2+2*M**2+ 3*s*v**2)/(s-M**2) - (4*m**2 + 2*m**2*(4*M**2+3*s*v**2+s ) + (M**2+s*v**2)**2 )/(s*v**2*(s-M**2))*np.log(1 + s*v**2/M**2) )

    return fc*( sigma_s + sigma_t + sigma_st)

def sigma_eff(v, g,M,m):
    v = v/c
    return (sigma_a(2*v, g, M, m) + sigma_r(2*v, g, M, m))/2


#### Total cross section, numerical #####
# t channel
def sigma_t_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (t**2 + 2*s*t + 2*s**2*v**2 + 8*m**4)/(t-M**2)**2
    return fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# u channel
def sigma_u_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (t**2 - 8*m**2*(s + t) + s**2 + 24*m**4)/(t + s*v**2 + M**2)**2
    return fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# t-u interference
def sigma_tu_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (4*m**4 - 2*(v**2 + 3)*m**2*s + s**2)/((t-M**2)*(t + s*v**2 + M**2))
    return - fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# s-channel
def sigma_s_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  (8*((2*m**2 + s*v**2 + t)**2 - s*(s-2*m**2) + (2*m**2 - t)**2) - 2*s*(4*(s-2*m**2) - 8*s))/(s - M**2)**2
    return fc *Prefactor * quad(integrand, -s*v**2, 0)[0]

# s-t interference
def sigma_st_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  - 8*(4*m**4 + 6*m**2*s*v**2 + 2*m**2*s + 8*m**2*t + s**2*v**4 + 2*s*v**2*t + t**2)/((s-M**2)*(t-M**2))
    return fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# Attractive
def sigma_a_Num(v, g, M, m):
    return sigma_t_num(v, g, M, m) + sigma_s_num(v, g, M, m) + sigma_st_num(v, g, M, m)

# Repulsive
def sigma_r_Num(v, g, M, m):
    return sigma_t_num(v, g, M, m) + sigma_u_num(v, g, M, m) + 2*sigma_tu_num(v, g, M, m)


# Total effective
def sigma_eff_Num(v, g, M, m):
    return (sigma_a_Num(v, g, M, m))# + sigma_r_Num(v, g, M, m) )/2

#### Transfer cross section, numerical #####
# t channel
def sigmaT_t_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(t**2 + 2*s*t + 2*s**2*v**2 + 8*m**4)/(t-M**2)**2
    return -2/(s*v**2) * fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# u channel
def sigmaT_u_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(t**2 - 8*m**2*(s + t) + s**2 + 24*m**4)/(t + s*v**2 + M**2)**2
    return -2/(s*v**2) * fc*Prefactor*quad(integrand, -s*v**2, 0, )[0]

# t-u interference
def sigmaT_tu_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = -(g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(4*m**4 - 2*(v**2 + 3)*m**2*s + s**2)/((t-M**2)*(t + s*v**2 + M**2))
    return  -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# s-channel
def sigmaT_s_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  t*(8*((2*m**2 + s*v**2 + t)**2 - s*(s-2*m**2) + (2*m**2 - t)**2) - 2*s*(4*(s-2*m**2) - 8*s))/(s - M**2)**2
    return -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]

# s-t interference
def sigmaT_st_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  - t*8*(4*m**4 + 6*m**2*s*v**2 + 2*m**2*s + 8*m**2*t + s**2*v**4 + 2*s*v**2*t + t**2)/((s-M**2)*(t-M**2))
    return -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# Attractive Transfer
def sigmaT_a_Num(v, g, M, m):
    return sigmaT_t_num(v, g, M, m) + sigmaT_s_num(v, g, M, m) + sigmaT_st_num(v, g, M, m)

# Repulsive Transfer
def sigmaT_r_Num(v, g, M, m):
    return sigmaT_t_num(v, g, M, m) + sigmaT_u_num(v, g, M, m) + 2*sigmaT_tu_num(v, g, M, m)


# Transfer effective
def sigmaT_Num(v, g, M, m):
    return (sigmaT_a_Num(v, g, M, m) + sigmaT_r_Num(v, g, M, m) )/2



#### Analytical cross section #####

#  Particle-Particle: chi chi -> chi chi
def sigma_r(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = (2*g**4/(d*np.pi*s*m))
    Term1 = ((2*s + 3*M**2)*s*v**2 + 2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2 + s*v**2))
    LogTerm1 = -(s*v**2*(3*M**2 + 4*m**2) + 2*(M**2 + 2*m**2)**2 - 4*m**4)/(s*v**2*(2*M**2 + s*v**2))
    return  fc*Prefactor*(Term1 + LogTerm1*np.log(1 + s*v**2/M**2) )


# Particle-Antiparticle: chi chibar -> chi chibar
def sigma_a(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = g**4/(d*np.pi*s*m)
    sigma_t = Prefactor*( (s*v**2*(2*s + 3*M**2)+2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2+s*v**2)) - (M**2 + s)/(s*v**2)*np.log(1 + s*v**2/M**2) )
    sigma_s = Prefactor/3*( (12*m**4 + 6*(v**2 + 1)*s*m**2 + s**2*v**4 ) / (s - M**2)**2 )
    sigma_st = -Prefactor/2*( (16*m**2+2*M**2+ 3*s*v**2)/(s-M**2) - (4*m**2 + 2*m**2*(4*M**2+3*s*v**2+s ) + (M**2+s*v**2)**2 )/(s*v**2*(s-M**2))*np.log(1 + s*v**2/M**2) )

    return fc*( sigma_s + sigma_t + sigma_st)

def sigma_eff(v, g,M,m):
    v = v/c
    return (sigma_a(2*v, g, M, m) + sigma_r(2*v, g, M, m))/2


# Ensemble average cross section
def sigv_integrand(v, v0, g, M, m):
    return sigma_eff(v, g, M, m)*v*np.exp(-0.5*v**2/v0**2)*v**2


def sigv_T(v0, g, M, m):
    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)

    Prefactor = 4*np.pi/((2*np.pi*v0**2)**1.5 * m)
    Integral = quad(sigv_integrand, 0.1, vmax, args=(v0, g, M, m))[0]
    return Prefactor*Integral
