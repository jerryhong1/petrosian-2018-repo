#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# blazars from the 3LAC catalog, with gamma data from 3FGL catalog and some optical data from SDSS DR6.
# in particular, we analyze the ~400 fsrq's out of the 3LAC catalog.


# In[0]: read and parse data; 999 sources with known redshift
from __future__ import division
import csv
import numpy as np
import quasars as quas
from quasars import QuasarData, Band
import matplotlib.pyplot as plt
import math
import scipy
import scipy.interpolate as interp
import numpy.polynomial.polynomial as poly
from astropy.modeling import models, fitting
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 14)
plt.rcParams['figure.figsize'] = (8.0, 6.0)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = 'true'
plt.rcParams['ytick.right'] = 'true'
plt.rcParams['errorbar.capsize'] =  2


names = [] # 3FGL names
cnames = [] #companion name
Z = [] # redshifts
Fr = [] # radio fluxes (mJy)
Fx = [] # x-ray fluxes (erg/s/cm^2)
RadTag = [] # radio flux source
V_USNO = [] # optical V magnitudes
V_SDSS = [] # optical V magnitudes
I_missing = [] # optical I magnitudes to be filled in

file_begin = 60
with open('3lac_high.tsv') as tsv:
    i = 0
    for line in csv.reader(tsv, delimiter="\t"):
        if (i >= file_begin):
            if 'fsrq' in line[7]:
                line = [s.replace(" ", "") for s in line]
                names.append(line[0])
                cnames.append(line[2])
                Z.append(line[8])
                Fr.append(line[9])
                RadTag.append(line[10])
                Fx.append(line[11])
                V_USNO.append(line[12])
                V_SDSS.append(line[13])
                I_missing.append(0.0)
#                find magnitudes on NED for objects without optical data provided
#                if(not line[12].strip() and not line[13].strip()): 
#                    print line[2]
        i = i + 1

# pull data from the larger catalog
names_3fgl = []
Fg_3fgl = []
Gamma_3fgl = []
file_begin = 79
with open('3fgl.tsv') as tsv:
    i = 0
    for line in csv.reader(tsv, delimiter="\t"):
        if (i >= file_begin):
            line = [s.replace(" ", "") for s in line]
            names_3fgl.append(line[2])
            Fg_3fgl.append(line[13])
            Gamma_3fgl.append(line[14])
        i = i + 1

Fg = [] # 
Gamma = [] # spectral indeces
for i in range(len(names)):
    index = [j for j in range(len(names_3fgl)) if names[i] in names_3fgl[j]][0]
    Fg.append(Fg_3fgl[index])
    Gamma.append(Gamma_3fgl[index])
    
# pull out SDSS i-band optical data (for ALL sources, not just those that are missing it)
file = open("3lac_high_opt_all.txt","r")
lines = file.readlines()
file.close()

for line in lines:
    linesplit = line.strip().split('|')
    linesplit = [s.replace(" ", "") for s in linesplit]
    name = linesplit[0]
    
    try:
        index = cnames.index(name)
    except ValueError:
        index = -1
    
    if(index != -1):
        if(linesplit[13]):
            I_missing[index] = float(linesplit[13])
            
# In[1]: identify where data is present/missing; 0.0 indicates missing

def str2float(Z):
    Z_ = []
    for z in Z:
        if(z):
            Z_.append(float(z))
        else:
            Z_.append(0.0)  
    return np.array(Z_)

Z = str2float(Z)
Fr = str2float(Fr)
Fx = str2float(Fx)
V_USNO = str2float(V_USNO)
V_SDSS = str2float(V_SDSS)
Fg = str2float(Fg) 
Fg = Fg * 1e-12 #over 100 MeV–100 GeV, erg cm−2 s−1, from power-law fit, 1 decimal place
Fr = Fr * 1e-26 #convert from mJy to cgs
Gamma = str2float(Gamma)
I_missing = np.array(I_missing)
#z_index = np.nonzero(Z)[0] #indeces where Z is known; will create a quasar object based on it

# In[]: radio band: doesn't really matter, except to correlate with gamma
fmin_r = 2.5e-26 #NVSS data only
k_radio = lambda z: quas.kcorrect(z, quas.alpha_rad)
r_band = Band('r', fmin_r, Fr, k_radio)

# In[3]: Gamma band: for Lmin, Zmax assuming average value for photon spectral index. 
fmin_g = 3e-12 # see 3LAC paper
Gamma_g = 2.44
alpha_g = 1 - Gamma_g
k_gamma = lambda z: quas.kcorrect(z, alpha_g)

def g_lum(z, f, gamma):
    alpha_g = gamma - 1    
    kc = (1 + Z[i])**(1 - alpha_g)
    return 4 * math.pi * (quas.d_lum(z)**2) * f / kc

#L = np.array([g_lum(Z[i], Fg[i], Gamma[i]) for i in range(len(Gamma))])
Fg = np.array([f if f > fmin_g else 0. for f in Fg ])
g_band = Band('\gamma', fmin_g, Fg, k_gamma)
#g_band.set_luminosity(L)
fsrq = QuasarData(Z, [g_band, r_band])
fsrq.sort()

# truncate everything below L = 1e45
trunc_index = np.where(g_band.L < 1e45)[0]
truncate = lambda A: np.array([A[i] if 1e45 < g_band.L[i] else 0.0 for i in range(len(A))])
#g_band.set_zmax(truncate(g_band.Zmax))
#g_band.set_luminosity(truncate(g_band.L))
#g_band.set_min_luminosity(truncate(g_band.Lmin))


# In[3]: Gamma band: defining a different Lmin for each source
L = np.array([g_lum(Z[i], Fg[i], Gamma[i]) for i in range(len(Gamma))])
Lmin = np.array([g_lum(Z[i], fmin_g, Gamma[i]) for i in range(len(Gamma))])
Lmin_avg = np.array([g_lum(Z[i], fmin_g, 2.44) for i in range(len(Gamma))])

#index = [i for i in range(len(fsrq.Z)) if Z[i] != 0 and L[i] != 0]
plt.figure()
plt.semilogy(Z, L, '.', markersize = 2)
plt.semilogy(Z, Lmin_avg, '.')
plt.semilogy(Z, Lmin, '.')
plt.title("Lgamma vs. z")
plt.savefig('lgamma.png')

plt.figure()
plt.errorbar(Z, np.log10(L), yerr = [np.log10(L) - np.log10(Lmin), np.zeros(len(Z))],
             fmt = 'o', markersize = 3, linewidth = 0.5)
plt.plot(fsrq.Z, np.log10(g_band.Lmin))

plt.figure()
plt.plot(np.log10(Lmin), np.log10(L), '.')
plt.plot(range(40, 50), range(40, 50))
plt.xlabel('Lmin')
plt.ylabel('L')


# In[]: The monotonicity problem

w = [i for i in range(len(fsrq.Z)) if Z[i] != 0 and L[i] != 0]
Zl = Z[w]
Ll, Lminl = quas.localize(Z[w], L[w], Lmin[w], 5.7)
i = 250

j1 = [m for m in range(0, len(Zl)) if (Ll[m] > Lminl[i] and Zl[m] < Zl[i])] #based off of L-z plot
j2 = [m for m in range(0, len(Zl)) if (Ll[m] > Lminl[i] and Lminl[m] < Lminl[i])] #based off of L-Lmin plot

plt.figure()
ax = plt.gca()
ax.set_yscale("log", nonposy='clip')
plt.title('Associated set based off of L-z plot')
plt.errorbar(Zl, Ll, yerr = [Ll - Lminl, np.zeros(len(Zl))],
             fmt = 'o', markersize = 3, linewidth = 0.5, color = 'black')
plt.errorbar(Zl[j1], Ll[j1], yerr = [Ll[j1] - Lminl[j1], np.zeros(len(Zl[j1]))],
             fmt = 'o', markersize = 3, linewidth = 0.5, color = 'blue')
plt.semilogy([0, Zl[i]], [Lminl[i], Lminl[i]], linewidth = 1, color = 'red')
plt.semilogy([Zl[i], Zl[i]], [Lminl[i] ,1e47], linewidth = 1, color = 'red')
plt.plot(Zl[i], Ll[i], '.', markersize = 12, color = 'red')

plt.figure()
plt.title('Associated set based off of L-Lmin plot')
plt.loglog(Lminl, Ll, '.', color = 'black')
plt.loglog(Lminl[j2], Ll[j2], '.', color = 'blue')
plt.semilogy([0, Lminl[i]], [Lminl[i], Lminl[i]], linewidth = 1, color = 'red')
plt.semilogy([Lminl[i], Lminl[i]], [Lminl[i] ,1e48], linewidth = 1, color = 'red')
plt.plot(np.linspace(10**43, 10**48, 10), np.linspace(10**43, 10**48, 10))

plt.figure()
ax = plt.gca()
ax.set_yscale("log", nonposy='clip')
plt.title('Associated set based off of L-Lmin plot')
plt.errorbar(Zl, Ll, yerr = [Ll - Lminl, np.zeros(len(Zl))],
             fmt = 'o', markersize = 3, linewidth = 0.5, color = 'black')
plt.errorbar(Zl[j2], Ll[j2], yerr = [Ll[j2] - Lminl[j2], np.zeros(len(Zl[j2]))],
             fmt = 'o', markersize = 3, linewidth = 0.5, color = 'blue')

# In[4]: Optical band: converting to flux
# 1. create separate bands with separate k-corrections (USNO)
# 2. create QuasarData fsrq_opt with both bands? 
# 3. merge all data to create new band (MUST ENSURE fsrq_opt.Z  = fsrq.Z)

def magtoflux(V, f0):
    F = []
    for v in V:
        if v == 0.0:
            F.append(0.0)
        else:
            F.append(quas.magtoflux(v, f0))
    return F

F_USNO = magtoflux(V_USNO, quas.f0_v)
F_SDSS = magtoflux(V_SDSS, quas.f0_v)
F_SDSS_NED = magtoflux(I_missing, quas.f0_i)

# assume V band is 550 nm: convert everything to 2500 A
F_USNO = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_USNO]
F_SDSS = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS]
F_SDSS_NED = [quas.bandtoband(f, quas.lambda_i, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS_NED]

fmin_USNO = quas.bandtoband(quas.magtoflux(21, quas.f0_v), quas.lambda_v, quas.lambda_opt, quas.alpha_opt)
fmin_SDSS = quas.bandtoband(0.08317e-26, quas.lambda_i, quas.lambda_opt, quas.alpha_opt)

# merge 3LAC and NED SDSS data (should add 6 data points): prioritize 3LAC SDSS data
F = []
for i in range(len(F_SDSS)):
    if F_SDSS[i] != 0.0:
        F.append(F_SDSS[i])
    else: F.append(F_SDSS_NED[i])

F_SDSS = F

# assume average k correction for USNO data:
k_opt_USNO = lambda z: quas.kcorrect(z, quas.alpha_opt)

o_USNO = Band('USNO', fmin_USNO, F_USNO, k_opt_USNO)
o_SDSS = Band('SDSS', fmin_SDSS, F_SDSS, quas.kcorrect_opt)

fsrq_o = QuasarData(Z, [o_USNO, o_SDSS])
fsrq_o.sort()

# In[]: truncate and merge data: we should only care about Zmax, L, Lmin. If SDSS and USNO data both present, favor SDSS.
plt.figure()
plt.semilogy(fsrq_o.Z, o_USNO.L, '.', markersize = 2, color = 'black')
plt.semilogy(fsrq_o.Z, o_USNO.Lmin, color = 'black')
plt.semilogy(fsrq_o.Z, o_SDSS.L, '.', markersize = 2, color = 'red')
plt.semilogy(fsrq_o.Z, o_SDSS.Lmin, color = 'red')

def merge(attr):
    a = []
    for i in range(len(fsrq_o.Z)):
        if o_SDSS.__dict__['L'][i] != 0.0 :
            a.append(o_SDSS.__dict__[attr][i])
        else: a.append(o_USNO.__dict__[attr][i])
    return np.array(a)

Zmax = merge('Zmax')
L = merge('L')
Lmin = merge('Lmin')

# truncate:
trunc_index = np.where(Lmin > L)[0]
o_band_trunc = Band('o', None, None, None)
o_band_trunc.set_zmax(Zmax[trunc_index])
o_band_trunc.set_luminosity(L[trunc_index])
o_band_trunc.set_min_luminosity(Lmin[trunc_index])


truncate = lambda A: np.array([A[i] if Lmin[i] < L[i] else 0.0 for i in range(len(A))])
Zmax = truncate(Zmax)
L = truncate(L)

plt.semilogy(fsrq_o.Z, L, '.', color = 'blue')
plt.semilogy(fsrq_o.Z, Lmin, color = 'blue')

o_band = Band('o', None, None, None)
o_band.set_zmax(Zmax)
o_band.set_luminosity(L)
o_band.set_min_luminosity(Lmin)

fsrq.addband(o_band, True) #true = already sorted
fsrq.addband(o_SDSS, True)
fsrq.addband(o_USNO, True)
fsrq_trunc = QuasarData(fsrq_o.Z[trunc_index], [o_band_trunc])


###################################
# In[]: correlation analysis
###################################

Alpha, R = fsrq.correlation_analysis(o_band, g_band)
# k_opt ≈ 3.5, k_gamma ≈ 5.5

r_band.set_k(4.8)
g_band.set_k(4.0)
Alpha_rg, R_rg = fsrq.rvsalpha(r_band, g_band)

# In[5]: Plots
figurepath = '../figures/fermi3lac-fsrq/'

# L_opt vs. z
plt.figure(figsize=(8,6))
ax = plt.gca()
ax.set_yscale("log", nonposy='clip')
index = [i for i in range(len(fsrq.Z)) if o_band.L[i] != 0 and g_band.L[i] != 0]
plt.errorbar(fsrq.Z[index], o_band.L[index],
             yerr = [o_band.L[index] - o_band.Lmin[index], np.zeros(len(fsrq.Z[index]))],
             fmt = 'o', markersize = 3, elinewidth = 0.3, color = 'black')
j = [i for i in range(len(fsrq_trunc.Z)) if o_band_trunc.L[i]]
plt.errorbar(fsrq_trunc.Z[j], o_band_trunc.L[j], capsize = 1,
             yerr = [np.zeros(len(fsrq_trunc.Z[j])), o_band_trunc.Lmin[j] - o_band_trunc.L[j]],
             fmt = 'o', markersize = 2, elinewidth = 0.3, color = 'red')
plt.title(r"$L_{opt}$ vs. $z$")

# L_gam vs. z
plt.figure(figsize=(8,6))
plt.semilogy(fsrq.Z[index], g_band.L[index], '.', markersize = 3, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.semilogy(fsrq.Z[nonindex], g_band.L[nonindex], '.', markersize = 2, color = 'red')
plt.semilogy(fsrq.Z[index], g_band.Lmin[index])
plt.title(r"$L_{\gamma}$ vs. $z$")

# L_opt vs. z, just limited set
plt.figure(figsize=(8,6))
ax = plt.gca()
ax.set_yscale("log", nonposy='clip')
index = o_band.limited_indeces
plt.errorbar(fsrq.Z[index], o_band.L[index], 
             yerr = [o_band.L[index] - o_band.Lmin[index], np.zeros(len(fsrq.Z[index]))],
             fmt = 'o', markersize = 3, linewidth = 0.5, color = 'black')
plt.title(r"$L_{opt}$ vs. $z$ (optically-limited set only)")

# L_gam vs. z, just limited set
plt.figure(figsize=(8,6))
index = g_band.limited_indeces
plt.semilogy(fsrq.Z[index], g_band.L[index], '.', markersize = 2, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.semilogy(fsrq.Z[index], g_band.Lmin[index])
plt.title(r"$L_{\gamma}$ vs. $z$ (gamma-limited set only)")


# tau vs k
sigma = [-1, 1]
for b in [fsrq.bands[0]]: #just gamma
    Taus = b.tau_array
    Ks = b.k_array
    k_g = b.k_g

    plt.figure()
    plt.title(r"$\tau$ vs $k_{" + b.name + "}$")
    plt.plot(Ks, Taus)
    axes = plt.gca()
    axes.set_xlim([min(Ks), max(Ks)])
    axes.set_ylim([-3,3])
    plt.xlabel(r"$k_{" + b.name + "}$", fontsize = 14)
    plt.ylabel(r"$\tau$", fontsize = 14)
    plt.minorticks_on()         

    plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
    plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
    plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)

    plt.plot([k_g, k_g], [-4, 0], color = 'red', linewidth = 1)
    plt.text(k_g + 0.01, 0.1, r"$k_{" + b.name + "} = $" + str(round(k_g, 2)), color = 'red', fontsize = 14)
    k_raderr = [np.interp(i, -Taus, Ks) for i in sigma]
    plt.plot([k_raderr[0], k_raderr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
    plt.plot([k_raderr[1], k_raderr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)
    
    plt.savefig(figurepath + '/tauvsk-fsrq-' + b.name + '.eps')
    plt.show()
    
# In[]: Plots con't: alpha vs r
o_band.set_k(3.5)
Alpha, R = fsrq.rvsalpha(o_band, g_band)

plt.figure()
plt.plot(Alpha, R, color = 'black', linewidth = 1)
#plt.plot(Alpha[::10], R[::10], marker = '.', markersize = 8)
alpha = np.interp(0, -np.array(R), Alpha)
plt.plot([alpha, alpha], [-4, 0], color = 'red', linewidth = 1)
plt.ylabel(r"Correlation Coefficient")
plt.xlabel(r"$\alpha$")
plt.title(r"Correlation of $L_{cr}^{\prime\; \gamma} = L_{\gamma}'(L_0/L_{opt}')^\alpha$ vs. $L_{opt}'$")
plt.text(alpha + 0.02, 0.01, r"$\alpha = $ " + str(round(alpha, 2)), color = 'red', fontsize = 14)

axes = plt.gca()
plt.plot(np.arange(-1, 3), np.zeros(4), color = 'black', linewidth = 1)
axes.set_xlim([min(Alpha), 1])
axes.set_ylim([min(R), max(R)])
plt.minorticks_on()
plt.savefig(figurepath + 'r-alpha-og.eps')


Alpha, R = Alpha_rg, R_rg
plt.figure()
plt.plot(Alpha, R, color = 'black', linewidth = 1)
#plt.plot(Alpha[::10], R[::10], marker = '.', markersize = 8)
alpha = np.interp(0, -np.array(R), Alpha)
plt.plot([alpha, alpha], [-4, 0], color = 'red', linewidth = 1)
plt.ylabel(r"Correlation Coefficient")
plt.xlabel(r"$\alpha$")
plt.title(r"Correlation of $L_{cr}^{\prime\; \gamma} = L_{\gamma}'(L_0/L_{rad}')^\alpha$ vs. $L_{rad}'$")
plt.text(alpha + 0.02, 0.01, r"$\alpha = $ " + str(round(alpha, 2)), color = 'red', fontsize = 14)

axes = plt.gca()
plt.plot(np.arange(-1, 3), np.zeros(4), color = 'black', linewidth = 1)
axes.set_xlim([min(Alpha), 1])
axes.set_ylim([min(R), max(R)])
plt.minorticks_on()
plt.savefig(figurepath + 'r-alpha-rg.eps')





###############################################################################
# In[4]: gamma density, luminosity functions
###############################################################################
b = g_band
i = b.limited_indeces
#i = [i for i in range(len(b.L))] 
L_raw = b.L[i] # raw, non-local luminosity
Lmin_raw = b.Lmin[i]
L, Lmin = quas.localize(fsrq.Z[i], b.L[i], b.Lmin[i], b.k_g)
#L, Lmin = quas.localize(fsrq.Z[i], b.L[i], b.Lmin[i], 4.0)
Z = fsrq.Z[i]
Zmax = g_band.Zmax[i]

plt.figure()
plt.semilogy(fsrq.Z[i], L, '.')
plt.semilogy(fsrq.Z[i], Lmin)

# In[]: Cumulative density function (sigma)
#Z_cdf = np.hstack((np.arange(0, 0.5, 0.02), np.arange(0.5, 3, 0.05)))
Z_cdf = Z
CDF = np.array([quas.cdf(z, Z, L, Lmin) for z in Z_cdf])
N = np.array([quas.N_z(z, Z, L) for z in Z_cdf])
CDF_raw = np.array([quas.cdf(z, Z, L_raw, Lmin_raw) for z in Z_cdf])

# In[]: Attempts to fit CDF

# CDF log-log fit
incr = [j for j in range(0, len(Z_cdf) - 10) if Z_cdf[j] != Z_cdf[j + 1]] #start at 1 for optical
#incr.append(-2)
Zfit = Z_cdf[incr]
CDFfit = CDF[incr]
Nfit = N[incr]

#w = np.hstack((np.zeros(5) + 30., np.zeros(20) + 0., np.zeros(20) + 5, np.zeros(len(Zfit) - 45) + 20)) # first 10 at weight 0.5 for optical

'''
f = interp.splrep(np.log10(1 + Zfit), np.log10(CDFfit), w, k = 3, s = 10)
fit = lambda x: interp.splev(x, f)
dfit = lambda x: interp.splev(x, f, der = 1)
'''

'''
p = poly.polyfit(np.log10(1 + Zfit), np.log10(CDFfit), 7, w = w)
fit = lambda x: poly.polyval(x, p)
dfit = lambda x: poly.polyval(x, poly.polyder(p))
'''

'''
def broken_power(x, A, xc, g1, g2):
    return A * x**g1 / (1 + (x / xc)**g2)
m = models.SmoothlyBrokenPowerLaw1D(200, 1, 1, 1, 1)
f = fitting.LevMarLSQFitter()
fit = f(m, Zfit, CDFfit)
dfit = lambda x: scipy.misc.derivative(fit, x, dx = 0.005)
'''

# regular fit:

fit = interp.UnivariateSpline(Zfit, CDFfit, k = 3, s = 10000)
dfit = fit.derivative()


fitN = interp.UnivariateSpline(Zfit, Nfit, k = 3, s = 1000)
dfitN = fitN.derivative()


'''
def broken_power(x, A, xc, g1, g2):
    return A * x**g1 / (1 + (x / xc)**g2)
p, foo = scipy.optimize.curve_fit(broken_power, Zfit, CDFfit)#, sigma = w**-1)
print p
fit = lambda x: broken_power(x, *p)
dfit = lambda x: scipy.misc.derivative(fit, x, dx = 0.005)
'''

# In[]: Plot CDF, fit, rho

logZcurve = np.arange(min(np.log10(1 + Zfit)), max(np.log10(1 + Zfit)), 0.02)
Zcurve = 10**logZcurve - 1

# CDF log plot
plt.figure()
plt.loglog(Z_cdf, CDF, '.',)
plt.loglog(Z_cdf, N, '.', Z_cdf, CDF_raw)
#plt.loglog(Zcurve, 10**fit(logZcurve)) #log log fit
plt.title(r"Cumulative Density Function $\log(\sigma(z))$ for $L_{" + b.name + '}$')
plt.xlabel(r"$z$")
#plt.xlabel(r"$Z = 1 + z$")
plt.ylabel(r"$\sigma$")
plt.minorticks_on()
plt.savefig(figurepath + '/CDF-N-' + b.name + '.eps')

# non log-log CDF plot
plt.figure()
plt.plot(Z_cdf, CDF, '.',  markersize = 2)
plt.plot(Z_cdf, N, '.',  markersize = 2)
#plt.plot(Zcurve, 10**fit(logZcurve)) # log-log fit
plt.plot(Zcurve, fit(Zcurve)) # regular fit
plt.plot(Zcurve, fitN(Zcurve))
plt.xlabel(r"$z$")
plt.ylabel(r"$\sigma(z)$")
plt.minorticks_on()

def loglogdevolution(z):
    Z = 1 + z
    logZ = np.log10(Z)
    sigma = 10**fit(logZ)
    return dfit(logZ) * sigma / Z / quas.dV_dz(z)

#dsigma_dz = [loglogdevolution(z) * quas.dV_dz(z) for z in Zfit]    
dsigma_dz = [dfit(z) for z in Zfit] #regular fit
dN_dz = [dfitN(z) for z in Zfit]
#rho = [loglogdevolution(z) for z in Zfit] #log-log fit
rho = [dfit(z) / quas.dV_dz(z) for z in Zfit] #regular fit
rho = np.array(rho) * 3.086e+24**3 #convert from /cm^3 to /Mpc^3

plt.figure()
plt.plot(Zfit, dsigma_dz, '.')
plt.plot(Zfit, dN_dz, '.')
plt.title(r"$ \frac{d \sigma}{dz}$ and $\frac{d N}{dz}$")
plt.xlabel(r"$z$")
plt.minorticks_on()
plt.savefig(figurepath + '/ds-dz-' + b.name + '.eps')

plt.figure()
plt.plot(Zfit, rho, '.')
plt.title(r"Density Evolution $\rho(z)$ for $L_{" + b.name + "}$'")
plt.xlabel(r"$z$")
plt.ylabel(r"$\rho(z)$")
plt.minorticks_on()
plt.savefig(figurepath + '/rho-' + b.name + '.eps')

'''
LDF = np.array([quas.ldf(Zfit[z], L, b.k_g, rho[z]) for z in range(len(Zfit))])
plt.figure()
plt.plot(Zfit, LDF, '.')
plt.title(r"Luminosity Density Function \textsterling$(z)$ for $L_{" + b.name + '}'$')
'''
# In[]: Cumulative Luminosity Function (Phi)
Z_zmax = np.hstack((np.arange(0.01, 4.99, 0.01), np.arange(5, 500, 0.5)))
log_Lmin = [np.log10(g_band.min_luminosity(z) / quas.g(z, g_band.k_g)) for z in Z_zmax]
x = lambda l: 0.0 if l == 0.0 else np.interp(np.log10(l), 
                                                     log_Lmin, Z_zmax, 0, float("inf"))
Zmax_local = np.array([x(l) for l in L])

CLF = np.array([quas.clf(l, fsrq.Z[i], L, Zmax_local) for l in L])
#    CLF= CLF / CLF[0] # NORMALIZE
N = np.array([quas.N_L(l, fsrq.Z[i], L) for l in L])
N_raw = np.array([quas.N_L(l, fsrq.Z[i], L_raw) for l in L_raw])
CLF_raw = np.array([quas.clf(l, fsrq.Z[i], L_raw, Zmax) for l in L_raw])

# In[]: Phi graph and LF calculations
plt.figure()
plt.loglog(L, CLF, '.', L, N, '.')
plt.title(r"Cumulative Luminosity Function $\Phi(L)$ and $N$ for $L_{" + b.name + "}'$")
plt.xlabel(r"$L_{" + b.name + '}$')
plt.ylabel(r"$\Phi(L_{" + b.name + '})$')
plt.minorticks_on()

Lfit = L[np.argsort(L)][2:]
CLFfit = CLF[np.argsort(L)][2:]
incr = [j for j in range(0, len(Lfit) - 1) if Lfit[j] != Lfit[j + 1] and Lfit[j] != 0] #eliminate repeat values of L
incr.append(-1)
Lfit = Lfit[incr]
CLFfit = CLFfit[incr]

fit = interp.UnivariateSpline(np.log10(Lfit), np.log10(CLFfit), s = .4)
dfit = fit.derivative()
Lcurve = np.linspace(min(Lfit), max(Lfit), 100)
plt.loglog(Lcurve, 10**fit(np.log10(Lcurve)))
plt.savefig(figurepath + '/CLF' + b.name + '.eps')

# raw (de-evolved) luminosity analysis
#plt.loglog(L_raw, CLF_raw, '.', L_raw, N_raw, '.')

LF = np.array([- 10**fit(np.log10(l)) * dfit(np.log10(l)) / l for l in Lfit])
plt.figure()
plt.loglog(Lfit[3:], LF[3:], '.')
plt.title(r"Luminosity Function $\psi(L)$ for $L_{" + b.name + '}$')
plt.xlabel(r"$L_{" + b.name + '}$')
plt.ylabel(r"$\psi(L_{" + b.name + '})$')
plt.savefig(figurepath + '/LF' + b.name + '.eps')
plt.minorticks_on()
    

# In[]: compare LF with Ajello and Singal:

L_singal = [1.2e45, 2e45, 3.1e45, 5.7e45, 8.3e45, 1.45e46, 2.4e46, 5e46]
LF_singal = [2e-7, 1e-7, 9e-8, 2.5e-8, 1.5e-8, 5e-9, 1e-9, 3e-10]
L_ajello = [6e45, 1e46, 2e46, 3e46, 6e46]
LF_ajello = [2e-8, 5e-9, 3e-9, 5e-10, 3e-10]    

plt.figure()
plt.loglog(L_singal, LF_singal, '.', L_ajello, LF_ajello, Lfit[7:-1], 10**36.5*LF[7:-1], '.')

# In[]: Various troubleshooting

plt.figure() 
plt.loglog((Z_zmax[0:500]), [quas.dV_dz(z) for z in Z_zmax[0:500]], '.')
plt.plot()
plt.title(r'$\frac{dV}{dz}$ vs. $z$')

plt.figure()
plt.semilogx(Fg, Gamma, '.')
#plt.gca().set_xlim([10-12, -8.7])
#plt.gca().set_ylim([1, 3.5])


#plt.figure()
#i = np.where((o_band.L != 0) & (g_band.L != 0))
#plt.plot(np.log10(o_band.L[i]), np.log10(g_band.L[i]), '.')
#plt.gca().set_xlim([0, 10])

plt.figure()
plt.plot(Zfit[2:], rho4[2:], 'b.', Zfit[2:], rho5[2:], 'r.')
plt.title(r'$\rho(z)$, blue: k = 4, red: k = 5.5')

plt.figure()
plt.loglog(Lfit[2:], LF4[2:], 'b.', Lfit[2:], LF5[2:], 'r.')
plt.title(r'$\psi(z)$, blue: k = 4, red: k = 5.5')