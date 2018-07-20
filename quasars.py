# coding: utf-8
#!/usr/bin/env python2
from __future__ import division
import numpy as np
from scipy.integrate import quad
import scipy.stats as stat
import math

class QuasarData:
    #Z = redshifts
    #bands = [Band]. assume they have been set up as defined in Band:
    #prerequisites: data has already been truncated, fluxes are fluxes of that band
    
    def __init__(self, Z, bands):
        print "initializing quasar data \n"
        self.Z = Z
        self.bands = bands
        
        for b in bands:
            print "band: " + b.name
            F = b.F
            
            L = np.array([b.luminosity(Z[i], F[i]) for i in range(len(Z))])
            b.set_luminosity(L)
            print "luminosity set up"
            
            Lmin = np.array([b.min_luminosity(Z[i]) for i in range(len(Z))])
            b.set_min_luminosity(Lmin)
            print "minimum luminosity set up"
            
            if(len(bands) > 1):
                Zmax = np.array([self.Zmax(L, b) for i in range(len(Z))])
                b.set_zmax(Zmax)
                print "ZMax set up"
            
            print "\n"
    
    def luminosity(self, z, f, band):
       return band.luminosity(z, f)
    
    def Zmax(self, L, band):
        Z = np.hstack((np.arange(0.001, 4.99, 0.001), np.arange(5, 500, 0.5)))
        log_Lmin = [np.log10(self.luminosity(z, band.fmin, band.k)) for z in Z]
        x = lambda logl: np.interp(logl, log_Lmin, Z, 0, float("inf"))
        return x(L)
    
    def size(self):
        return len(self.Z)
    
    #sort all data by redshift
    def sort(self):
        index = np.argsort(self.Z)
        self.Z = self.Z[index]
        
        for b in self.bands:
            L = b.L[index]
            b.set_luminosity(L)
            Lmin = b.Lmin[index]
            b.set_min_luminosity(Lmin)    
            if(len(self.bands) > 1):
                Zmax = b.Zmax[index]
                b.set_zmax(Zmax)    
                
    
class Band:
    fmin = 0.0
    k = lambda x: x
    
    def __init__(self, bandname, fmin, F, k_correction):
        self.name = bandname
        self.fmin = fmin
        self.F = F
        self.k = k_correction
        
    def luminosity(self, z, f):
        kc = self.k(z)
        return 4 * math.pi * (d_lum(z)**2) * f / kc

    def min_luminosity(self, z):
        kc = self.k(z)
        return 4 * math.pi * (d_lum(z)**2) * self.fmin / kc
    
    def set_luminosity(self, L):
        self.L = L
        
    def set_min_luminosity(self, Lmin):
        self.Lmin = Lmin
        
    def set_zmax(self, Zmax):
        self.Zmax = Zmax
        

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
#k-corrections for richard et al. (2006), which are normalized at z = 2
f = open('kcorr_quas.txt', 'r')
K = []
Z_k = []
for line in f:
    Z_k.append(line[0:4])
    K.append(line[5:11])
Z_k = [float(z) for z in Z_k]
K = [float(k) for k in K]

#print K, Z_k
alpha_opt = -0.5
def k_opt(z): #with sign convention: m_intrinsic = m_obs - K, and L = 4 pi d^2 f/k(z)
    k_avg = -2.5 * (1 + alpha_opt) * math.log10(1 + z) 
    if (z > max(Z_k)):
        k = k_avg  #assume no emission line effect
    else:  
        k = np.interp(z, Z_k, K) - 2.5 * (1 + alpha_opt) * math.log10(1 + 2)
    return 10**(-k / 2.5)


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
        if (len(j) == 1 or len(j) == 2): continue
        
        #determine rank
        L_rank = stat.rankdata(L_ass, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (len(L_ass))
        
        resid = resid + (rank - exp)
        var = var + (1/12.0 * (-1 + (len(L_ass)-1)**2))
        
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

#1-dimensinoal tau vs. k
def tauvsk(Z, L, Lmin, K):
    s = []
    for k in K:
        L_loc, Lmin_loc = localize(Z, L, Lmin, k)
        print('\nk = ' +  str(k))
        t = tau(Z, L_loc, Lmin_loc)
        print ('\ntau = ' + str(t))
        s.append(t)
    return np.array(s)
