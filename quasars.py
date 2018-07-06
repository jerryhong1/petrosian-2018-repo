#!/usr/bin/env python2
import numpy as np
from scipy.integrate import quad
import scipy.stats as stat
import math

def convert(Z, M, m_max, m_min, f0):
    F = [magtoflux(m, f0) for m in M]
    L = [lum(Z[i], F[i]) for i in range(0,len(F))]
    fmin = magtoflux(m_min)
    Lmin = [lmin(L[i], F[i], fmin) for i in range(0, len(F))]

    F = np.array(F)
    L = np.array(L)
    Lmin = np.array(Lmin)
        
    data = np.stack((Z, M, F, L, Lmin))
    data = np.transpose(data)
    data, data_trunc = truncate(data, m_max, m_min)
    return data, data_trunc    

#luminosity distance in cm. assume the standard ΛCDM cosmology: H0 = 71 km s−1 Mpc−1, ΩΛ = 0.7, and Ωm = 0.3
#work in cgs
def integrand(x, a, b):
    return 1 / math.sqrt(a + b*(1+x)**3)

def d_lum(z):
    Omega_l = 0.7
    Omega_m = 0.3
    cH0 = 1.30291e+28 #cm, ratio of c/H0. H0 = 71 km/s/Mpc
    integral = quad(integrand, 0, z, args = (Omega_l, Omega_m))
    return integral[0] * (1+z) * cH0

#Convert magnitudes to fluxes (in cgs). see Abazajian et al. (2009) §3.3
# (multiply Jy by e-23 to get cgs: ergs/s/cm^2/Hz)
def magtoflux(m,f0):
#    f0 = 3.631e-20 value for Singal (2013)
    return f0*10**(-m/2.5)


#calculate luminosity, minimum luminosity
#k-corrections for richard et al. (2006)
f = open('kcorr_quas.txt', 'r')
K = []
Z_k = []
for line in f:
    Z_k.append(line[0:4])
    K.append(line[5:11])

def k(z): #with sign convention: m_intrinsic = m_obs - K, and L = 4 pi d^2 f/k(z)
    k = float(K[int(100 * z)])            
    return 10**(-k / 2.5)

#TODO: convert from one wavelength to another
def lum(z, f):
    kc = k(z)
    return 4*math.pi*(d_lum(z)**2)*f/kc

#fmin = 0.083e-26 #value for Singal et al. (2013)
def lmin(l, f, fmin):
    return l*fmin/f


#concatentate quasar data?
def concatenate(Z, M, F, L, Lmin):
    data = np.stack((Z, M, F, L, Lmin))
    data = np.transpose(data)
    return data

#double truncate data based off of apparent magnitude
def truncate(data, m_max, m_min):
    data_temp = [];
    data_trunc = [];
    for i in range(0, len(data[:,1])):
        if(data[i,1] > m_min):
            if(data[i,1] < m_max):
                data_temp.append(data[i,:])
            else:
                data_trunc.append(data[i,:])
    
    data_trunc = np.array(data_trunc)        
    data = np.array(data_temp)
    return data, data_trunc

#tau calculation
def tau(Z, L, Lmin): #as defined in singal, petrosian papers in 2010s. tau = (∑resid)/(sqrt(∑variance))
    resid = 0
    var = 0
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])] #see petrosian
        j.append(i)
        L_ass = L[j]
        
        #determine rank
        L_rank = stat.rankdata(L_ass, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (1 + len(L_ass))
        
        resid = resid + (rank - exp)
        var = var + (1/12.0 * (-1 + len(L_ass)**2))
        
        #troubleshooting
        if(i % 500 == 0): print i, resid, var
        
    t = resid / math.sqrt(var)
    return t

#Local Luminosity creation, given k:
#g(z) such that L' = L/g(z). Taken from Singal et al. (2013) and considers high-redshift objects.
def g(z, k):
    z_cr = 3.7
    return (1 + z)**k /(1 + ((1 + z) / z_cr)**k)

def localize(Z, L, Lmin, k):
    L_local = [L[p] / g(Z[p],k) for p in range(0, len(L))]
    Lmin_local = [Lmin[p] / g(Z[p],k) for p in range(0, len(L))]
    return np.array(L_local), np.array(Lmin_local)

def tauvsk(K, data):
    s = []
    for k in K:
        L_loc, Lmin_loc = localize(data[:,0], data[:,3], data[:,4], k)
        print('\nk = ' +  str(k))
        t = tau(data[:,0], L_loc, Lmin_loc)
        print ('\ntau = ' + str(t))
        s.append(t)
    return np.array(s)
