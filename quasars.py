# coding: utf-8
#!/usr/bin/env python2
from __future__ import division
import numpy as np
from scipy.integrate import quad
import scipy.stats as stat
import math

# this file defines:
# 1. QuasarData, a class to store all of a quasar sample's data and analyze it,
#    including each band, local luminosity determination, correlation between bands
# 2. Band, a class to store all the quasar data related to a band of light: fluxes, minimum flux
#    luminosities, minimum luminosities, k correction formula, sources which are "band-limited",
#    local luminosity parameters.
# 3. Commonly used formulas in analyzing quasar data.

class QuasarData:
    # Z = redshifts
    # bands = [Band]. assume they have been set up as defined in Band:
    # prerequisites: data has already been truncated, fluxes are fluxes of that band

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

            if not hasattr(b, 'ZMax'):
                if(len(bands) > 1):
                    Zmax = b.Zmax(b.L)
                    b.set_zmax(Zmax)
                    print "ZMax set up"

            print "\n"

    def luminosity(self, z, f, band):
       return band.luminosity(z, f)

    # returns number of quasars in sample
    def size(self):
        return len(self.Z)

    # sort all data by redshift: DO BEFORE doing correlation analysis
    # so that limited_indeces isn't messed up
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

#     for 2 bands, divides dataset into 2 disjoint "band-limited" sets. (accepts 2 bandnames as input)
#     prompts user for k values to try tau analysis.
#     prompts use for k array to determine k where tau = 0;
#     saves this value of k, as well as the k-tau array for graphing
#     finally, localizes luminosities and graphs alpha vs. rad
    def correlation_analysis(self, bandname1, bandname2):
        band1 = self.bands[np.where(bandname1 == np.array(self.bandnames))[0][0]]
        band2 = self.bands[np.where(bandname2 == np.array(self.bandnames))[0][0]]
        band1.add_index(band2, np.where(band1.Zmax < band2.Zmax)[0])
        band2.add_index(band1, np.where(band2.Zmax < band1.Zmax)[0])

        self.tauvsk(band1, band2)
        self.tauvsk(band2, band1)

        return self.rvsalpha(band1, band2)

    def tauvsk(self, band1, band2): #find k_g of band1 when analyzing w.r.t. band2
        b = band1
        index = b.limited_indeces[band2]
        print "beginning tau analysis, band " + b.name + ", " + str(len(index)) + " objects"

        k = raw_input("select a value of k to try, or press enter to skip: ")
        while(k):
            k = float(k)
            tauvsk(self.Z[index], b.L[index], b.Lmin[index], [k])
            k = raw_input("select a value of k to try, or press enter to skip: ")

        k = raw_input("select a range (\"min, max\") of k to analyze, or press enter to skip: ")    
        if(k):
            K = [float(i) for i in k.split(", ")]
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
        L_l1, foo = localize(self.Z, band1.L, band1.Lmin, band1.k_g)
        L_l2, foo = localize(self.Z, band2.L, band2.Lmin, band2.k_g)

        Alpha = np.arange(0,1,0.005)
        R = rvsalpha(L_l1, L_l2, Alpha)
        return Alpha, R


class Band:

    def __init__(self, bandname, fmin, F, k_correction):
        self.name = bandname
        self.fmin = fmin
        self.F = F
        self.kcorr = k_correction

        # dicitonary for indeces, assuming that there may be more than two bands which will be
        # all pairwise analyzed, in which case there will be different sets of "___-limited" for
        # each pairs of bands analyzed.
        self.limited_indeces = {}
        self.k_g = 0.0 #used for local luminosity analysis L/g(z)

    def luminosity(self, z, f):
        kc = self.kcorr(z)
        return 4 * math.pi * (d_lum(z)**2) * f / kc

    def min_luminosity(self, z):
        kc = self.kcorr(z)
        return 4 * math.pi * (d_lum(z)**2) * self.fmin / kc

    def Zmax(self, L):
        Z = np.hstack((np.arange(0.001, 4.99, 0.001), np.arange(5, 500, 0.5)))
        log_Lmin = [np.log10(self.min_luminosity(z)) for z in Z]
        x = lambda logl: np.interp(logl, log_Lmin, Z, 0, float("inf"))
        return x(np.log10(L))

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

    def add_index(self, band, index):
        self.limited_indeces[band] = index


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
#    f0 = 3.631e-20 value for i band, Singal (2013)
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
        var = var + (1/12.0 * (1 + (len(L_ass)-1)**2))

        #troubleshooting
        if(i % 500 == 0): print i, resid, var

    t = resid / math.sqrt(var)
    return t

#Local Luminosity creation, given k:
#g(z) such that L' = L/g(z). Taken from Singal et al. (2013); considers high-redshift objects.
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

def rvsalpha(L1, L2, Alpha):
    L0 = 1e30
    LumCr = lambda L1, L2, alpha: L2 * (L0/L1)**alpha
    R = []
    for a in Alpha:
        LCr = [LumCr(L1[i], L2[i], a) for i in range(len(L1))]
        r = stat.linregress(np.log10(L1), np.log10(LCr)).rvalue
        R.append(r)
    return R