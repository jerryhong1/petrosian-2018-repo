#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# blazars from the 3LAC catalog, with gamma radiation data from 3FGL catalog. 
# in particular, the ~400 fsrq's.


# In[0]: read and parse data; 999 sources with known redshift
from __future__ import division
import csv
import numpy as np
import quasars as quas
from quasars import QuasarData, Band
import matplotlib.pyplot as plt
import math
import scipy.interpolate as interp
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
            if ('fsrq' in line[7]):
                line = [s.replace(" ", "") for s in line]
                names.append(line[0])
                cnames.append(line[2])
                Z.append(line[8])
                Fr.append(line[9])
                Fx.append(line[11])
                V_USNO.append(line[12])
                V_SDSS.append(line[13])
                I_missing.append(0.0)
                #find magnitudes on NED for objects without optical data provided
#                if(line[12].strip() and line[13].strip()): 
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

Fg = []
Gamma = []
for i in range(len(names)):
    index = [j for j in range(len(names_3fgl)) if names[i] in names_3fgl[j]][0]
    Fg.append(Fg_3fgl[index])
    Gamma.append(Gamma_3fgl[index])
    
# pull out missing optical data: take i band
file = open("3lac_high_opt.txt","r")
lines = file.readlines()
file.close()

for line in lines:
    linesplit = line.strip().split('|')
    linesplit = [s.replace(" ", "") for s in linesplit]
    name = linesplit[1]
    
    try:
        index = cnames.index(name)
    except ValueError:
        index = -1
    
    if(index != -1):
        if(linesplit[16]):
            I_missing[index] = float(linesplit[16])
            
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
#z_index = np.nonzero(Z)[0] #indeces where Z is known; will create a quasar object based on it

# In[3]: Gamma band: assuming average value for photon spectral index

def k_g(z):    
    Gamma_g = 2.45
    alpha_g = Gamma_g - 1
    return (1 + z)**(1 - alpha_g)

fmin_g = 2e-12 # see 3LAC paper
g_band = Band('g', fmin_g, Fg, k_g)
fsrq = QuasarData(Z, [g_band])
fsrq.sort()

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

f0 = 3.631e-20  #assuming singal numbers again

# merge SDSS and USNO optical data
V = []
for i in range(len(Z)):
    if(V_USNO[i] != 0.0):
        V.append(V_USNO[i])
    else:
        V.append(V_SDSS[i])

Fo = []
for v in V:
    if v == 0.0:
        Fo.append(0.0)
    else:
        Fo.append(quas.magtoflux(v, f0))
Fo = np.array(Fo)


# add missing optical data
for i in range(len(I_missing)):
    if I_missing[i] != 0.0:
        f = quas.magtoflux(I_missing[i], f0)
        Fo[i] = quas.bandtoband(f, quas.lambda_i, quas.lambda_v, quas.alpha_opt)
        
        
Fo = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in Fo]

fmin_o = 0.08317e-26 * (quas.lambda_opt / 7625)**(-quas.alpha_opt)
#fmin_o = 0.02317e-26 * (quas.lambda_opt / 7625)**(-quas.alpha_opt)
#fmin_o = 0

#truncate data
Fo_ = []
Fo_trunc = []
Z_trunc = []
for i in range(len(Fo)):
    if Fo[i] > fmin_o:
        Fo_.append(Fo[i])
    else:
        Fo_.append(0.0)
        Z_trunc.append(Z[i])
        Fo_trunc.append(Fo[i])
        
Fo = Fo_

o_band = Band('o', fmin_o, Fo, quas.k_opt)
o_band_trunc = Band('o', fmin_o, Fo_trunc, quas.k_opt)
fsrq.addband(o_band)
fsrq_trunc = QuasarData(Z_trunc, [o_band_trunc])

# In[]: correlation analysis

Alpha, R = fsrq.correlation_analysis('o', 'g')


# In[3.5]: weird shit going on with tau; probably has to do with the g(z) when k is large

Ll, Lminl = quas.localize(fsrq.Z, g_band.L, g_band.Lmin, g_band.k_g)
plt.figure()
plt.plot(fsrq.Z, np.log10(Ll), '.', markersize = 2)
plt.plot(fsrq.Z, np.log10(Lminl))
plt.title("L/g(z, k) vs. z where tau = 0")
plt.savefig('local.png')

# In[3]: Gamma band: defining a different Lmin for each source

def lum(z, f, gamma):
    alpha_g = gamma - 1
    
    kc = (1 + Z[i])**(1 - alpha_g)
    return 4 * math.pi * (quas.d_lum(z)**2) * f / kc

L = np.array([lum(Z[i], Fg[i], Gamma[i]) for i in range(len(Gamma))])
Lmin = np.array([lum(Z[i], 2e-12, Gamma[i]) for i in range(len(Gamma))])

index = [i for i in range(len(fsrq.Z)) if Z[i] != 0 and L[i] != 0]
plt.figure()
plt.plot(Z[index], np.log10(L[index]), '.', markersize = 2)
plt.plot(fsrq.Z[index], np.log10(g_band.Lmin[index]))
plt.plot(Z[index], np.log10(Lmin[index]), '.')
plt.title("Lgamma vs. z")
plt.savefig('lgamma.png')

plt.figure()
plt.errorbar(Z[index], np.log10(L[index]), yerr = [np.log10(L[index]) - np.log10(Lmin[index]), np.zeros(len(Z))], fmt = 'o', markersize = 3, linewidth = 0.5)

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

# In[5]: Plots
figurepath = '../figures/fermi3lac-fsrq/'

#L_opt vs. z
plt.figure()
index = [i for i in range(len(fsrq.Z)) if o_band.L[i] != 0]
plt.plot(fsrq.Z[index], np.log10(o_band.L[index]), '.', markersize = 2, color = 'black')
plt.plot(fsrq.Z[index], np.log10(o_band.Lmin[index]))
plt.plot(fsrq_trunc.Z, np.log10(o_band_trunc.L), '.', markersize = 1, color = 'red')
plt.title(r"$L_{opt}$ vs. $z$")

#L_gam vs. z
plt.figure()
plt.plot(fsrq.Z[index], np.log10(g_band.L[index]), '.', markersize = 2, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.plot(fsrq.Z[nonindex], np.log10(g_band.L[nonindex]), '.', markersize = 1, color = 'red')
plt.plot(fsrq.Z, np.log10(g_band.Lmin))
plt.title(r"$L_{\gamma}$ vs. $z$")

#L_opt vs. z, just limited set
plt.figure()
index = o_band.limited_indeces
plt.plot(fsrq.Z[index], np.log10(o_band.L[index]), '.', markersize = 2, color = 'black')
plt.plot(fsrq.Z[index], np.log10(o_band.Lmin[index]))
plt.title(r"$L_{opt}$ vs. $z$ (optically-limited set only)")

#L_gam vs. z, just limited set
plt.figure()
index = g_band.limited_indeces
plt.plot(fsrq.Z[index], np.log10(g_band.L[index]), '.', markersize = 2, color = 'black')
nonindex = [i for i in range(len(fsrq.Z)) if o_band.L[i] == 0]
plt.plot(fsrq.Z, np.log10(g_band.Lmin))
plt.title(r"$L_{\gamma}$ vs. $z$ (gamma-limited set only)")



#tau vs k
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

#alpha vs r
plt.figure()
plt.plot(Alpha, R, color = 'black', linewidth = 1)
#plt.plot(Alpha[::10], R[::10], marker = '.', markersize = 8)
alpha = np.interp(0, -np.array(R), Alpha)
plt.plot([alpha, alpha], [-4, 0], color = 'red', linewidth = 1)
plt.ylabel(r"Correlation Coefficient")
plt.xlabel(r"$\alpha$")
plt.title(r"Correlation of $L_{cr}^{\prime\; rad} = L_{rad}'(L_0/L_{opt}')^\alpha$ vs. $L_{opt}'$")
plt.text(alpha + 0.02, 0.01, r"$\alpha = $ " + str(round(alpha, 2)), color = 'red', fontsize = 14)

axes = plt.gca()
plt.plot(np.arange(-1, 3), np.zeros(4), color = 'black', linewidth = 1)
axes.set_xlim([min(Alpha), 1])
axes.set_ylim([min(R), max(R)])
plt.minorticks_on()
plt.savefig(figurepath + 'r-alpha.eps')

# In[4]: density, luminosity functions

for b in fsrq.bands:
    Z_cdf = np.arange(0.05, 3.2, 0.05)
    i = b.limited_indeces
    L, Lmin = quas.localize(fsrq.Z[i], b.L[i], b.Lmin[i], b.k_g)
    
#    plt.figure()
#    plt.plot(fsrq.Z[i], np.log10(L), '.')
#    plt.plot(fsrq.Z[i], np.log10(Lmin))
    
    CDF = np.array([quas.cdf(z, fsrq.Z[i], L, Lmin) for z in Z_cdf])
#    fit = interp.InterpolatedUnivariateSpline(Z_cdf, CDF, k = 3)
    
    plt.figure()
    plt.plot(Z_cdf, CDF, '.')
#    plt.plot(Z_cdf, fit(Z_cdf))
    plt.title(r"Cumulative Density Function $\sigma(z)$ for $L_{" + b.name + '}$')
    
    rho = np.array([quas.devolution(z, Z_cdf, CDF) for z in Z_cdf])
    plt.figure()
    plt.plot(Z_cdf, rho)
    plt.title(r"Density Evolution $\rho(z)$ for $L_{" + b.name + '}$')
    
    LDF = np.array([quas.ldf(Z_cdf[z], L, b.k_g, rho[z]) for z in range(len(Z_cdf))])
    plt.figure()
    plt.plot(Z_cdf, LDF)
    plt.title(r"Luminosity Density Function \textsterling$(z)$ for $L_{" + b.name + '}$')
    
    L_lf = np.arange(min(np.log10(L)), max(np.log10(L)), 0.05)
    L_lf = 10**L_lf
    CLF = np.array([quas.clf(l, fsrq.Z[i], L, Lmin) for l in L_lf])
    plt.figure()
    plt.plot(np.log10(L_lf), np.log10(CLF), '.')
    plt.title(r"Cumulative Luminosity Function $\Phi(L)$ for $L_{" + b.name + '}$')
    
    LF = np.array([quas.lf(l, L_lf, CLF) for l in L_lf])
    plt.figure()
    plt.plot(np.log10(L_lf), np.log10(LF), '.')
    plt.title(r"Luminosity Function $\psi(L)$ for $L_{" + b.name + '}$')

#plt.figure()
#i = np.where((o_band.L != 0) & (g_band.L != 0))
#plt.plot(np.log10(o_band.L[i]), np.log10(g_band.L[i]), '.')
#plt.gca().set_xlim([0, 10])