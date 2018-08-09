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
import scipy.interpolate as interp
import numpy.polynomial.polynomial as poly
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 14)
plt.rcParams['figure.figsize'] = (8.0, 6.0)

names = [] # 3FGL names
cnames = [] #companion name
Z = [] # redshifts
Fr = [] # radio fluxes (mJy)
Fx = [] # x-ray fluxes (erg/s/cm^2)
V_USNO = [] # optical V magnitudes
V_SDSS = [] # optical V magnitudes
I_missing = [] # optical I magnitudes to be filled in

file_begin = 60
with open('3lac_high.tsv') as tsv:
    i = 0
    for line in csv.reader(tsv, delimiter="\t"):
        if (i >= file_begin):
            if 'fsrq' in line[7]: # and not (quas.isFloat(line[8]) and float(line[8]) < 0.03): # take out first data point
                line = [s.replace(" ", "") for s in line]
                names.append(line[0])
                cnames.append(line[2])
                Z.append(line[8])
                Fr.append(line[9])
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
Gamma = str2float(Gamma)
I_missing = np.array(I_missing)
#z_index = np.nonzero(Z)[0] #indeces where Z is known; will create a quasar object based on it

# In[3]: Gamma band: for Lmin, Zmax assuming average value for photon spectral index. 

def k_g(z):    
    Gamma_g = 2.45
    alpha_g = Gamma_g - 1
    return (1 + z)**(1 - alpha_g)

def g_lum(z, f, gamma):
    alpha_g = gamma - 1
    
    kc = (1 + Z[i])**(1 - alpha_g)
    return 4 * math.pi * (quas.d_lum(z)**2) * f / kc

#L = np.array([g_lum(Z[i], Fg[i], Gamma[i]) for i in range(len(Gamma))])
fmin_g = 3e-12 # see 3LAC paper
Fg = np.array([f if f > fmin_g else 0 for f in Fg ])
g_band = Band('\gamma', fmin_g, Fg, k_g)
#g_band.set_luminosity(L)
fsrq = QuasarData(Z, [g_band])
fsrq.sort()

# truncate everything below L = 1e45
#trunc_index = np.where(g_band.L < 1e45)[0]
#truncate = lambda A: np.array([A[i] if 1e45 < g_band.L[i] else 0.0 for i in range(len(A))])
#g_band.set_zmax(truncate(g_band.Zmax))
#g_band.set_luminosity(truncate(g_band.L))
#g_band.set_min_luminosity(truncate(g_band.Lmin))


'''
K =  np.arange(5,6,0.1)
Tau = quas.tauvsk(fsrq.Z[index], g_band.L[index], g_band.Lmin[index], K)
k_gam = np.interp(0, -Tau, K)
plt.figure()
plt.plot(K, Tau)
plt.title("Tau vs k, assuming avg Gamma")
plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)
axes = plt.gca()
axes.set_xlim([min(K), max(K)])
axes.set_ylim([-3,3])
plt.savefig('tauvsk.png')
'''

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

# assume V band is 550 nm: convert everything to optical
F_USNO = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_USNO]
F_SDSS = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS]
F_SDSS_NED = [quas.bandtoband(f, quas.lambda_i, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS_NED]

fmin_USNO = quas.bandtoband(quas.magtoflux(21, quas.f0_v), quas.lambda_v, quas.lambda_opt, quas.alpha_opt)
fmin_SDSS = quas.bandtoband(0.08317e-26, quas.lambda_i, quas.lambda_opt, quas.alpha_opt)

# merge 3LAC and NED SDSS data (should add 6 data points): prioritize...3LAC SDSS data? (ask petrosian)
F = []
for i in range(len(F_SDSS)):
    if F_SDSS[i] != 0.0:
        F.append(F_SDSS[i])
    else: F.append(F_SDSS_NED[i])

F_SDSS = F

# assume average k correction for USNO data:
k_opt_USNO = lambda z: (1 + z)**(1 + quas.alpha_opt)

o_USNO = Band('USNO', fmin_USNO, F_USNO, k_opt_USNO)
o_SDSS = Band('SDSS', fmin_SDSS, F_SDSS, quas.k_opt)

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
trunc_index = np.where(Lmin[i] > L[i])[0]
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

# In[3]: Gamma band: defining a different Lmin for each source
fmin_g = 3e-12
L = np.array([g_lum(Z[i], Fg[i], Gamma[i]) for i in range(len(Gamma))])
Lmin = np.array([g_lum(Z[i], fmin_g, Gamma[i]) for i in range(len(Gamma))])
Lmin_avg = np.array([g_lum(Z[i], fmin_g, 2.45) for i in range(len(Gamma))])

#index = [i for i in range(len(fsrq.Z)) if Z[i] != 0 and L[i] != 0]
plt.figure()
plt.semilogy(Z, L, '.', markersize = 2)
plt.semilogy(Z, Lmin_avg, '.')
plt.semilogy(Z, Lmin, '.')
plt.title("Lgamma vs. z")
plt.savefig('lgamma.png')

plt.figure()
plt.errorbar(Z, np.log10(L), yerr = [np.log10(L) - np.log10(Lmin), np.zeros(len(Z))], fmt = 'o', markersize = 3, linewidth = 0.5)
plt.plot(fsrq.Z, np.log10(g_band.Lmin))

plt.figure()
plt.plot(np.log10(Lmin), np.log10(L), '.')
plt.plot(range(40, 50), range(40, 50))
plt.xlabel('Lmin')
plt.ylabel('L')

'''
Tau = quas.tauvsk(fsrq.Z[index], g_band.L[index], g_band.Lmin[index], K)
plt.figure()
plt.plot(K, Tau)
plt.title("Tau vs k, accounting for Gamma")
plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)
axes = plt.gca()
axes.set_xlim([min(K), max(K)])
axes.set_ylim([-3,3])
#plt.savefig('tauvsk.png')
'''

# In[]: correlation analysis

Alpha, R = fsrq.correlation_analysis(o_band, g_band)
# k_opt ≈ 3.5, k_gamma ≈ 5.5

# In[5]: Plots
figurepath = '../figures/fermi3lac-fsrq/'

# L_opt vs. z
plt.figure(figsize=(8,6))
index = [i for i in range(len(fsrq.Z)) if o_band.L[i] != 0 or g_band.L[i] != 0]
plt.semilogy(fsrq.Z[index], o_band.L[index], '.', markersize = 2, color = 'black')
plt.semilogy(fsrq.Z[index], o_band.Lmin[index])
plt.semilogy(fsrq_trunc.Z, o_band_trunc.L, '.', markersize = 1, color = 'red')
plt.title(r"$L_{opt}$ vs. $z$")

# L_gam vs. z
plt.figure(figsize=(8,6))
plt.semilogy(fsrq.Z[index], g_band.L[index], '.', markersize = 2, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.semilogy(fsrq.Z[nonindex], g_band.L[nonindex], '.', markersize = 1, color = 'red')
plt.semilogy(fsrq.Z[index], g_band.Lmin[index])
plt.title(r"$L_{\gamma}$ vs. $z$")

# L_opt vs. z, just limited set
plt.figure(figsize=(8,6))
index = o_band.limited_indeces
plt.semilogy(fsrq.Z[index], o_band.L[index], '.', markersize = 2, color = 'black')
plt.semilogy(fsrq.Z[index], o_band.Lmin[index])
plt.title(r"$L_{opt}$ vs. $z$ (optically-limited set only)")

# L_gam vs. z, just limited set
plt.figure(figsize=(8,6))
index = g_band.limited_indeces
plt.semilogy(fsrq.Z[index], g_band.L[index], '.', markersize = 2, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.semilogy(fsrq.Z, g_band.Lmin)
plt.title(r"$L_{\gamma}$ vs. $z$ (gamma-limited set only)")


# tau vs k
sigma = [-1, 1]
for b in fsrq.bands:
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
plt.savefig(figurepath + 'r-alpha.eps')

# In[4]: gamma density, luminosity functions

b = g_band
i = range(len(b.L))
L_raw = b.L[i] # raw, non-local luminosity
Lmin_raw = b.Lmin[i]
L, Lmin = quas.localize(fsrq.Z[i], b.L[i], b.Lmin[i], b.k_g)
#    L, Lmin = quas.localize(fsrq.Z[i], b.L[i], b.Lmin[i], 5.5)
Z = fsrq.Z[i]
Zmax = g_band.Zmax[i]

plt.figure()
plt.semilogy(fsrq.Z[i], L, '.')
plt.semilogy(fsrq.Z[i], Lmin)

# In[]: Cumulative density function (sigma)
CDF = np.array([quas.cdf(z, Z, L, Lmin) for z in Z])
N = np.array([quas.N_z(z, Z, L) for z in Z])
CDF_raw = np.array([quas.cdf(z, Z, L_raw, Lmin_raw) for z in Z])

# In[]: Plot CDF, determine Density Evolution (rho)
plt.figure()
plt.loglog(1 + Z, CDF, '.')
#plt.loglog(1 + Z, N, '.', 1 + Z, CDF_raw)
plt.title(r"Cumulative Density Function $\log(\sigma(z))$ for $L_{" + b.name + '}$')
plt.xlabel(r"$Z = 1 + z$")
plt.ylabel(r"$\sigma$")
plt.minorticks_on()

# CDF log-log fit
incr = [j for j in range(0, len(Z) - 1) if Z[j] != Z[j + 1]] #start at 1 for optical
Zfit = Z[incr]
CDFfit = CDF[incr]
w = np.hstack((np.zeros(4) + 10, np.zeros(1) + 0.5, np.arange(len(Zfit) - 5) + 1)) # first 10 at weight 0.5 for optical, 30 for gamma lol


fit = interp.UnivariateSpline(np.log10(1 + Zfit), np.log10(CDFfit), w, k = 3, s = 500)
dfit = fit.derivative()

'''
p = poly.polyfit(np.log10(1 + Zfit), np.log10(CDFfit), 9)
fit = lambda x: poly.polyval(x, p)
dfit = lambda x: poly.polyval(x, poly.polyder(p))
'''

'''
# regular fit:
fit = interp.UnivariateSpline(Zfit, CDFfit, s = 10000)
dfit = fit.derivative()
'''

logZcurve = np.arange(min(np.log10(1 + Zfit)), max(np.log10(1 + Zfit)), 0.02)
Zcurve = 10**logZcurve - 1
plt.loglog(1 + Zcurve, 10**fit(logZcurve)) #log log fit
plt.minorticks_on()


#non log-log fit check
plt.figure()
plt.plot(Z, CDF, '.')
plt.plot(Zcurve, 10**fit(logZcurve)) # log-log fit
#plt.plot(Zcurve, fit(Zcurve)) # regular fit
plt.xlabel(r"$z$")
plt.ylabel(r"$\sigma(z)$")
plt.minorticks_on()

def devolution(z):
    Z = 1 + z
    logZ = np.log10(Z)
    sigma = 10**fit(logZ)
    return dfit(logZ) * sigma / Z / quas.dV_dz(z)
    
rho = [devolution(z) for z in Zfit] #log-log fit
#rho = [dfit(z) / quas.dV_dz(z) for z in Zfit] #regular fit
plt.figure()
plt.plot(Zfit[1:], rho[1:], '.')
plt.title(r"Density Evolution $\rho(z)$ for $L_{" + b.name + "}$'")
plt.xlabel(r"$z$")
plt.ylabel(r"$\rho(z)$")
plt.minorticks_on()

'''
LDF = np.array([quas.ldf(Zfit[z], L, b.k_g, rho[z]) for z in range(len(Zfit))])
plt.figure()
plt.plot(Zfit, LDF, '.')
plt.title(r"Luminosity Density Function \textsterling$(z)$ for $L_{" + b.name + '}'$')
'''

# In[]: Cumulative Luminosity Function (Phi) and LF (psi)
CLF = np.array([quas.clf(l, fsrq.Z[i], L, Zmax) for l in L])
#    CLF= CLF / CLF[0] # NORMALIZE
N = np.array([quas.N_L(l, fsrq.Z[i], L) for l in L])
plt.figure()
plt.loglog(L, CLF, '.', L, N, '.')
plt.title(r"Cumulative Luminosity Function $\Phi(L)$ and $N$ for $L_{" + b.name + "}'$")
plt.xlabel(r"$L_{" + b.name + '}$')
plt.ylabel(r"$\Phi(L_{" + b.name + '})$')
plt.minorticks_on()

Lfit = L[np.argsort(L)][2:]
CLFfit = CLF[np.argsort(L)][2:]
incr = [j for j in range(1, len(Lfit) - 1) if Lfit[j] != Lfit[j + 1]] #eliminate repeat values of L
Lfit = Lfit[incr]
CLFfit = CLFfit[incr]

fit = interp.UnivariateSpline(np.log10(Lfit), np.log10(CLFfit), s = 10)
dfit = fit.derivative()
Lcurve = np.linspace(min(Lfit), max(Lfit), 100)
plt.loglog(Lcurve, 10**fit(np.log10(Lcurve)))

# raw (de-evolved) luminosity analysis
N_raw = np.array([quas.N_L(l, fsrq.Z[i], L_raw) for l in L_raw])
CLF_raw = np.array([quas.clf(l, fsrq.Z[i], L_raw, Zmax) for l in L_raw])
plt.loglog(L_raw, CLF_raw, '.', L_raw, N_raw, '.')

LF = np.array([- 10**fit(np.log10(l)) * dfit(np.log10(l)) / l for l in Lfit])
plt.figure()
plt.loglog(Lfit[3:-1], LF[3:-1], '.')
plt.title(r"Luminosity Function $\psi(L)$ for $L_{" + b.name + '}$')
plt.xlabel(r"$L_{" + b.name + '}$')
plt.ylabel(r"$\psi(L_{" + b.name + '})$')
plt.minorticks_on()

'''
Nfit = N[np.argsort(L)][2:]
fit = interp.UnivariateSpline(Lfit, Nfit, s = 1000)
plt.loglog(Lfit, Nfit, '.', Lfit, fit(Lfit))
plt.figure()
plt.loglog(Lfit, -fit.derivative()(Lfit), '.')
'''
    

# In[]: compare LF with Ajello and Singal:

L_singal = [1.2e45, 2e45, 3.1e45, 5.7e45, 8.3e45, 1.45e46, 2.4e46, 5e46]
LF_singal = [2e-7, 1e-7, 9e-8, 2.5e-8, 1.5e-8, 5e-9, 1e-9, 3e-10]
L_ajello = [6e45, 1e46, 2e46, 3e46, 6e46]
LF_ajello = [2e-8, 5e-9, 3e-9, 5e-10, 3e-10]    

plt.figure()
plt.loglog(L_singal, LF_singal, '.', L_ajello, LF_ajello, Lfit[7:-1], 10**37*LF[7:-1], '.')

 # In[]: Various troubleshooting

plt.figure() 
plt.plot(Zfit, [quas.dV_dz(z) for z in Zfit], '.')
plt.plot()

plt.figure()
plt.semilogx(Fg, Gamma, '.')
#plt.gca().set_xlim([10-12, -8.7])
#plt.gca().set_ylim([1, 3.5])


#plt.figure()
#i = np.where((o_band.L != 0) & (g_band.L != 0))
#plt.plot(np.log10(o_band.L[i]), np.log10(g_band.L[i]), '.')
#plt.gca().set_xlim([0, 10])
    