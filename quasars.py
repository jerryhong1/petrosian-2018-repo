# coding: utf-8
#!/usr/bin/env python2
from __future__ import division
import numpy as np
from scipy.integrate import quad
import scipy.stats as stat
import matplotlib.pyplot as plt
import math
import cPickle as pickle

# this file defines:
# 1. QuasarData, a class to store all of a quasar sample's data and analyze it,
#    including each band, local luminosity determination, correlation between bands
# 2. Band, a class to store all the quasar data related to a band of light: fluxes, minimum flux
#    luminosities, minimum luminosities, k correction formula, sources which are "band-limited",
#    local luminosity parameters.
# 3. Commonly used formulas in analyzing quasar data.

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

class QuasarData:
    # Z = redshifts
    # bands = [Band]. assume they have been set up as defined in Band. 
    # (Assume we are comparing two bands at once at most?)
    # prerequisites: data has already been truncated, fluxes are fluxes of that band, 
    # all bands are unique.
    # allows for missing data = 0.0

    # if bands have not had luminosity/Zmax calculated, calculates them based on the fluxes provided.
    def __init__(self, Z, bands):
        print "initializing quasar data \n"
        self.Z = Z
        self.bands = bands
        self.bandnames = [b.name for b in bands]


        for b in bands:
            print "band: " + b.name
            F = b.F

            if not hasattr(b, 'L'):
                L = np.array([b.luminosity(Z[i], F[i]) for i in range(len(Z))])
                b.set_luminosity(L)
                print "luminosity set up"

            if not hasattr(b, 'Lmin'):
                Lmin = np.array([b.min_luminosity(Z[i]) for i in range(len(Z))])
                b.set_min_luminosity(Lmin)
                print "minimum luminosity set up"

            if not hasattr(b, 'Zmax'):
                if(len(bands) <= 1):
                    a = raw_input("only 1 band. calculate Zmax? (not recommended if data is large) (y/n)" )
                    while (a != 'y' and a!= 'n'):
                        a = raw_input("only 1 band. calculate Zmax? (not recommended if data is large) (y/n)" )
                    if(a == 'y'):
                        Zmax = b.Zmax_calc(b.L)
                        b.set_zmax(Zmax)
                        print "Zmax set up"
                else:
                    Zmax = b.Zmax_calc(b.L)
                    b.set_zmax(Zmax)
                    print "Zmax set up"
                    
            print "\n"
            
    def addband(self, band): #or, overwrite.
        b = band
        Z = self.Z
        
        if (band.name not in self.bandnames):
            print "adding new band!"
            self.bands.append(b)
            self.bandnames.append(b.name)

        print "band: " + b.name
        
        if hasattr(self, 'sort_index'): #already sorted
            self.sortband(band) #should just sort b.F out, unless others are already calculated.
            
        F = b.F

        if not hasattr(b, 'L'):
            L = np.array([b.luminosity(Z[i], F[i]) for i in range(len(Z))])
            b.set_luminosity(L)
            print "luminosity set up"

        if not hasattr(b, 'Lmin'):
            Lmin = np.array([b.min_luminosity(Z[i]) for i in range(len(Z))])
            b.set_min_luminosity(Lmin)
            print "minimum luminosity set up"
        
        if not hasattr(b, 'Zmax'):
            Zmax = b.Zmax_calc(b.L)
            b.set_zmax(Zmax)
            print "Zmax set up"

        print "\n"
        
        
    def luminosity(self, z, f, band):
       return band.luminosity(z, f)

    # returns number of quasars in sample
    def size(self):
        return len(self.Z)

    # sort all data by redshift: DO BEFORE doing correlation analysis
    # so that limited_indeces isn't messed up
    def sort(self):
        if not hasattr(self, 'sort_index'):
            index = np.argsort(self.Z)
            self.sort_index = index # save index for later bands to be automatically sorted.
            self.Z = self.Z[index]
            print index

        for b in self.bands:
            self.sortband(b)
            
    def sortband(self, band):
        i = self.sort_index
        b = band
        
        if hasattr(b, 'L'):
            L = b.L[i]
            b.set_luminosity(L)
            
        if hasattr(b, 'Lmin'):
            Lmin = b.Lmin[i]
            b.set_min_luminosity(Lmin)
        
        if hasattr(b, 'Zmax'):
            Zmax = b.Zmax[i]
            b.set_zmax(Zmax)
        
        
#     for 2 bands, divides dataset into 2 disjoint "band-limited" sets. (accepts 2 bandnames as input)
#     prompts user for k values to try tau analysis.
#     prompts use for k array to determine k where tau = 0;
#     saves this value of k, as well as the k-tau array for graphing
#     finally, localizes luminosities and graphs alpha vs. rad
    def correlation_analysis(self, bandname1, bandname2):
        band1 = self.bands[np.where(bandname1 == np.array(self.bandnames))[0][0]]
        band2 = self.bands[np.where(bandname2 == np.array(self.bandnames))[0][0]]
        band1.add_index(np.where((band1.Zmax < band2.Zmax) & (band1.Zmax != 0))[0])
        band2.add_index(np.where((band2.Zmax < band1.Zmax) & (band2.Zmax != 0))[0])

        plt.figure()
        index = band1.limited_indeces
        plt.plot(self.Z[index], np.log10(band1.L[index]), '.', markersize = 2)
        plt.plot(self.Z[index], np.log10(band1.Lmin[index]))
        
        plt.figure()
        index = band2.limited_indeces
        plt.plot(self.Z[index], np.log10(band2.L[index]), '.', markersize = 2)
        plt.plot(self.Z[index], np.log10(band2.Lmin[index]))

    
        self.tauvsk(band1, band2)
        self.tauvsk(band2, band1)

        return self.rvsalpha(band1, band2)

    def tauvsk(self, band1, band2): #find k_g of band1 when analyzing w.r.t. band2
        b = band1
        index = b.limited_indeces
        print "beginning tau analysis, band " + b.name + ", " + str(len(index)) + " objects"

        k = raw_input("select a value of k to try, or press enter to skip: ")
        while(k):
            if (isFloat(k)):
                k = float(k)
                tauvsk(self.Z[index], b.L[index], b.Lmin[index], [k])
            k = raw_input("select a value of k to try, or press enter to skip: ")

        k = raw_input("select a range (\"min, max\") of k to analyze, or press enter to skip: ")    
        if(k):
            K = [float(i) for i in k.split(",")]
            while len(K) != 2:
                k = raw_input("select a range (\"min, max\") of k to analyze, or press enter to skip: ")
                K = [float(i) for i in k.split(", ")]

            K = np.arange(K[0], K[1], 0.1)

            Tau = tauvsk(self.Z[index], b.L[index], b.Lmin[index], K)

            b.set_k(np.interp(0, -Tau, K)) #assuming Tau is decreasing
            b.set_tau_array(Tau)
            b.set_k_array(K)
        
        print "k = " + str(b.k_g)

    def rvsalpha(self, band1, band2):
        i = np.where((band1.L != 0) & (band2.L != 0))[0]
        print i
        L_l1, foo = localize(self.Z[i], band1.L[i], band1.Lmin[i], band1.k_g)
        L_l2, foo = localize(self.Z[i], band2.L[i], band2.Lmin[i], band2.k_g)
        
        plt.figure()
        plt.plot(np.log10(L_l1), np.log10(L_l2), '.', markersize = 2)
        
        Alpha = np.arange(0,1,0.005)
        R = rvsalpha(L_l1, L_l2, Alpha)
        return Alpha, R


class Band:
    def __init__(self, bandname, fmin, F, k_correction):
        self.name = bandname
        self.fmin = fmin
        self.F = F
        self.kcorr = k_correction
#        self.L = []
#        self.Lmin = []

        # dicitonary for indeces, assume just two bands are being compared.
        self.limited_indeces = []
        self.k_g = 0.0 #used for local luminosity analysis L/g(z)

    def luminosity(self, z, f):
        kc = self.kcorr(z)
        return 4 * math.pi * (d_lum(z)**2) * f / kc

    def min_luminosity(self, z):
        kc = self.kcorr(z)
        return 4 * math.pi * (d_lum(z)**2) * self.fmin / kc

    def Zmax_calc(self, L):
        Z = np.hstack((np.arange(0.001, 4.99, 0.001), np.arange(5, 500, 0.5)))
        log_Lmin = [np.log10(self.min_luminosity(z)) for z in Z]
        x = lambda l: 0.0 if l == 0.0 else np.interp(np.log10(l), 
                                                     log_Lmin, Z, 0, float("inf"))
        return np.array([x(l) for l in L])

    # can accept empty values as 0.0
    def set_luminosity(self, L):
        self.L = L

    def set_min_luminosity(self, Lmin):
        self.Lmin = Lmin

    def set_zmax(self, Zmax):
        self.Zmax = Zmax

    def set_k(self, k):
        self.k_g = k

    def set_k_array(self, K):
        self.k_array = K

    def set_tau_array(self, Tau):
        self.tau_array = Tau

    def add_index(self, index):
        self.limited_indeces = index


################################################################################
#### MISCELLANEOUS, IMPORTANT FUNCTIONS
################################################################################
    
# luminosity distance in cm. assume the standard ΛCDM cosmology: 
# H0 = 71 km s−1 Mpc−1, ΩΛ = 0.7, and Ωm = 0.3
# work in cgs
def integrand(x, a, b):
    return 1 / math.sqrt(a + b*(1+x)**3)

Omega_l = 0.7
Omega_m = 0.3
cH0 = 1.30291e+28 #cm, ratio of c/H0. H0 = 71 km/s/Mpc
def r_comoving(z):
    integral = quad(integrand, 0, z, args = (Omega_l, Omega_m))
    return integral[0] * cH0

def d_lum(z):
    return r_comoving(z) * (1 + z)

# Convert magnitudes to fluxes (in cgs). see Abazajian et al. (2009) §3.3
# (multiply Jy by e-23 to get cgs: ergs/s/cm^2/Hz)
def magtoflux(m,f0):
#    f0 = 3.631e-20 value for i band, Singal (2013)
    return f0*10**(-m/2.5)


#Band to band conversion (SDSS bands): all wavelengths reported in angstroms
alpha_opt = -0.5    
alpha_opt = -0.6
lambda_v = 5510. 
lambda_r = 6231.
lambda_i = 7625
lambda_opt = 2500.


# alpha defined as F \propto \nu^\alpha 
# from band 1 to band 2
def bandtoband(f, lambda1, lambda2, alpha):
    return f *(lambda2 / lambda1)**(-alpha_opt)


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


# tau calculation
# as defined in singal, petrosian papers in 2010s. tau = (∑resid)/(sqrt(∑variance))
def tau(Z, L, Lmin): 
    resid = 0
    var = 0
    for i in range(len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])] # and Lmin[m] < Lmin[i])] #see petrosian
        j.append(i)
        size = len(j)
        L_ass = L[j]

        if (size <= 1): continue
    
        #determine rank
        L_rank = stat.rankdata(L_ass, 'max') #ranks of all luminosities
        rank = L_rank[-1] #determine rank of data point i, including data point in set.
        exp = (size + 1) / 2.

        resid = resid + (rank - exp)
        var = var + (-1 + size**2) / 12.0

        #troubleshooting
        if(i % 500 == 0): print i, resid, var

    t = resid / math.sqrt(var)
    return t

#Local Luminosity creation, given k:
#g(z) such that L' = L/g(z). 
# Taken from Singal et al. (2013); considers high-redshift objects.
def g(z, k):
    z_cr = 3.7
    return (1 + z)**k /(1 + ((1 + z) / z_cr)**k)

def g_(z, k):
    z_cr = 2.5
    return (1 + z)**k /(1 + ((1 + z) / z_cr)**k)


def localize(Z, L, Lmin, k):
    L_local = [L[p] / g(Z[p],k) for p in range(0, len(L))]
    Lmin_local = [Lmin[p] / g(Z[p],k) for p in range(0, len(L))]
    return np.array(L_local), np.array(Lmin_local)

#tau vs. k
def tauvsk(Z, L, Lmin, K):
    s = []
    for k in K:
        L_loc, Lmin_loc = localize(Z, L, Lmin, k)
        print('\nk = ' +  str(k))
        t = tau(Z, L_loc, Lmin_loc)
        print ('\ntau = ' + str(t))
        s.append(t)
    return np.array(s)

# after localizing luminosities, use correlation-reduced lumninosity 
# to calculate correlation
def rvsalpha(L1, L2, Alpha):
    L0 = 1e30
    LumCr = lambda L1, L2, alpha: L2 * (L0/L1)**alpha
    R = []
    for a in Alpha:
        LCr = [LumCr(L1[i], L2[i], a) for i in range(len(L1))]
        r = stat.linregress(np.log10(L1), np.log10(LCr)).rvalue
        R.append(r)
    return R

# cumulative density function:
def cdf(z, Z, L, Lmin):
    sigma = 1
    for i in range(len(Z)):
        #create associated sets
        if Z[i] > z: continue
        j = [m for m in range(0, len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])] # and Lmin[m] < Lmin[i])] #see petrosian
        size = len(j)
        if size > 0:
            sigma = sigma * (1. + 1. / size)
            
    return sigma
    
#rate of change of co-moving volume w.r.t. redshift
dV_dz = lambda z: 4 * math.pi * (r_comoving(z) ** 2) * cH0 * integrand(z, Omega_l, Omega_m)  # units of cm^3

# density evolution function: 
def devolution(z, Z, CDF):
    dsdz = np.interp(z, Z, np.gradient(CDF, Z))
    dVdz = dV_dz(z)
    return dsdz * 1. / dVdz

# energy production rate w.r.t. redshift
def ldf(z, L_local, k, rho):
    return sum(L_local) / len(L_local) * g(z, k) * rho * dV_dz(z)

# cumulative luminosity function
def clf(l, Z, L, Lmin):
    phi = 1
    for i in range(len(Z)):
        if L[i] < l: continue
        # note that luminosity must be greater than L[i], not Lmin[i]
        j = [m for m in range(0, len(Z)) if (L[m] > L[i] and Z[m] < Z[i])] 
        size = len(j)
        if size > 0:
            phi = phi * (1. + 1. / size)
    return phi

# luminosity function
def lf(l, L, CLF):
    return -np.interp(l, L, np.gradient(CLF, L))
    