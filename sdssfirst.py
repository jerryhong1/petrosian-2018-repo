#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import quasars as quas
from quasars import Band, QuasarData
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# In[1]: import file and put into QuasarData anLarad band objects.
file = open("SDSSFIRSTComplete.dat","r") 
lines = file.readlines()
file.close() 
data = np.empty((len(lines) - 5, 7))
for i in range(5, len(lines)):
    linearr = np.array(lines[i].strip().split())
    for j in range(0, len(linearr)):
        data[i - 5, j] = float(linearr[j])

# z = 0, d_lum = 1, k = 2, f_2500A = 3, f_1.4GHz = 4, L_o = 5, L_r = 6
#convert from mJy to cgs:
data[:,4] = data[:,4] * 1e-26

# In[2]: Determine minimum values
fmin_i = 0.083e-26 #see Singal et al. (2013)
alpha_opt = -0.5    
lambda_i = 7625. 
lambda_opt = 2500.
i_to_opt = lambda i: i *(lambda_opt / lambda_i)**(-alpha_opt)
fmin_opt = i_to_opt(fmin_i)

alpha_rad = -0.6
fmin_rad = 1.0e-26 #1 mJy
def k_rad(z):
    return (1 + z)**(1 + alpha_rad)

# In[3]: Create Objects
o_band = Band('o', fmin_opt, data[:,3], quas.k_opt)
o_band.set_luminosity(data[:,5])
r_band = Band('r', fmin_rad, data[:,4], k_rad)
r_band.set_luminosity(data[:,6])
sdssfirst = QuasarData(data[:,0], [o_band, r_band])

# In[4]: check plots
plt.plot(np.log10(sdssfirst.Z), np.log10(r_band.L), '.', markersize = 1)
plt.plot(np.log10(sdssfirst.Z), np.log10(r_band.Lmin), '-', color = 'red')
Z = np.hstack((np.arange(0.001, 4.99, 0.001), np.arange(5, 10, 0.5)))
plt.plot(np.log10(Z), np.log10([r_band.min_luminosity(z) for z in Z]), '-', color = 'red')

# In[5]: correlation analysis: tau vs k + r vs alpha
K_opt = np.arange(2.5, 4., 0.1)
K_rad = np.arange(4.5, 5.5, 0.1)
Alpha, R = sdssfirst.correlation_analysis(o_band, r_band)

# In[3]: PLOTs of log(L) vs. z
figurepath = '../figures/sdssxfirst'
plt.figure()
plt.loglog(data[:,5], data[:,6], '.', markersize=1, label="my data", color='black')
plt.xlabel(r"$\log(L_{opt})$")
plt.ylabel(r"$\log(L_{rad})$")
plt.title("Optical vs. Radio")
plt.minorticks_on()

#logl_opt vs z
plt.figure(figsize=(8,6))
plt.semilogy(data[:,0], data[:,5],'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10([LMinOpt(z) for z in data[:,0]]), linewidth = 1)
plt.xlabel(r"$z$", fontsize=14)
plt.ylabel(r"$L_{opt}$", fontsize=14)
plt.title(r"Luminosity at 2500 \r{A} vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([1e29, 1e33])
plt.minorticks_on()
plt.savefig(figurepath + 'OptSDSSFIRSTlogLz.png')

'''
#Test Zmax
i = 400
plt.plot(data[i,0], np.log10(data[i,5]), '.', markersize=12, color = 'red')
plt.plot(ZMaxOpt[i], np.log10(data[i,5]), '.', markersize=12, color = 'red')
'''

#logl_rad vs z
plt.figure(figsize=(8,6))
plt.plot(data[:,0], np.log10(data[:,6]),'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10([LMinRad(z) for z in data[:,0]]), linewidth = 1)
plt.xlabel(r"$z$", fontsize=14)
plt.ylabel(r"$L_{rad}$", fontsize=14)
plt.title(r"Luminosity at 1.4GHz vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
plt.minorticks_on()
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([29, 36.5])
plt.savefig(figurepath + 'RadSDSSFIRSTlogLz.png')
plt.show()

'''
#Test Zmax
plt.plot(data[i,0], np.log10(data[i,6]), '.', markersize=12, color = 'red')
plt.plot(ZMaxRad[i], np.log10(data[i,6]), '.', markersize=12, color = 'red')
'''
# In[] PLOT: Zmax Definition


plt.figure(figsize=(8,6))
#logl_rad vs z
plt.semilogy(data[:,0], r_band.L,'.', markersize=1, label="my data", color='#A0A0A0')
plt.semilogy(data[:,0], [r_band.min_luminosity(z) for z in data[:,0]], 
             linewidth = 1, color = 'red')
plt.xlabel("z", fontsize=14)
plt.ylabel(r"$L_{rad}$", fontsize=14)
plt.minorticks_on()
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([1e29, 3.16e36])

#Zmax
i = 2002
plt.semilogy(data[i,0], data[i,6], '.', markersize=12, color = 'black')
plt.semilogy([data[i,0],r_band.Zmax[i]], [r_band.L[i], r_band.L[i]],
             '--', color = 'black', linewidth = 2)
plt.semilogy([r_band.Zmax[i],r_band.Zmax[i]], [1, r_band.L[i]], '--', 
             color = 'black', linewidth = 2)
axes.text(data[i,0], data[i,6] * 10**0.3, r'\textbf{Source $\mathbf{i}$}',
        horizontalalignment='center',
        verticalalignment='center', 
        color = 'black', fontsize = 14,
        bbox=dict(facecolor='white', alpha=0.7, ec = 'none'))
axes.text(r_band.Zmax[i] + 0.1, 1e30, r'$z_{max, i}^{r}$',
        horizontalalignment='left',
        verticalalignment='top', 
        color = 'black', fontsize = 18)
plt.savefig(figurepath + 'zmax.png')
plt.show()

# In[] PLOT of tau vs k
sigma = [-1, 1]
#tau vs k_opt
Tau_opt = o_band.tau_array
K_opt = o_band.k_array
k_opt = o_band.k_g

plt.figure(3)
plt.title(r"$\tau$ vs $k_{opt}$ (SDSS X FIRST set)")
plt.plot(K_opt, Tau_opt)
axes = plt.gca()
axes.set_xlim([min(K_opt), max(K_opt)])
axes.set_ylim([-3,3])
plt.xlabel(r"$k_{opt}$", fontsize = 14)
plt.ylabel(r"$\tau$", fontsize = 14)

plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)

plt.plot([k_opt, k_opt], [-4, 0], color = 'red', linewidth = 1)
plt.text(k_opt + 0.01, 0.1, r"$k_{opt} = $ " + str(round(k_opt, 2)), 
         color = 'red', fontsize = 14)
k_opterr = [np.interp(i, -Tau_opt, K_opt) for i in sigma]
plt.plot([k_opterr[0], k_opterr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_opterr[1], k_opterr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)

plt.savefig(figurepath + 'tauvsk-first-opt.eps')
plt.show()

#tau vs k_rad
Tau_rad = r_band.tau_array
K_rad = r_band.k_array
k_rad = r_band.k_g

plt.figure(4)
plt.title(r"$\tau$ vs $k_{rad}$ (SDSS X FIRST set)")
plt.plot(K_rad, Tau_rad)
axes = plt.gca()
axes.set_xlim([min(K_rad), max(K_rad)])
axes.set_ylim([-3,3])
plt.xlabel(r"$k_{rad}$", fontsize = 14)
plt.ylabel(r"$\tau$", fontsize = 14)

plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)

plt.plot([k_rad, k_rad], [-4, 0], color = 'red', linewidth = 1)
plt.text(k_rad + 0.01, 0.1, r"$k_{rad} = $ " + str(round(k_rad, 2)), 
         color = 'red', fontsize = 14)
k_raderr = [np.interp(i, -Tau_rad, K_rad) for i in sigma]
plt.plot([k_raderr[0], k_raderr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_raderr[1], k_raderr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)
         
plt.savefig(figurepath + 'tauvsk-first-rad.eps')
plt.show()

# In[8]: More PLOTs

#Lopt' vs Lrad' (local)
plt.figure(figsize=(8,6))
L_optl, foo = quas.localize(sdssfirst.Z, o_band.L, o_band.Lmin, o_band.k_g)
L_radl, foo = quas.localize(sdssfirst.Z, r_band.L, r_band.Lmin, r_band.k_g)
plt.loglog(L_optl, L_radl, '.', markersize=1)
plt.xlabel(r"$L_{opt}'$", fontsize = 16)
plt.ylabel(r"$L_{rad}'$", fontsize = 16)
#plt.title("Local luminosity scatter plot")
axes = plt.gca()
plt.text(.74, 0.96, r"$g(z) = \frac{(1 + z)^k}{1 + (\frac{1 + z}{Z_{cr}})^k}$",
         horizontalalignment='left', verticalalignment='top',
         color = 'black', fontsize = 18, transform=axes.transAxes)
plt.text(.53, 0.96, r"$k_{rad} = $ " + str(round(k_rad, 2)) + "\n" + r"$k_{opt} = $ " 
         + str(round(k_opt, 2)) + "\n" + r"$Z_{cr} = $ " + str(3.7), color = 'black',
         horizontalalignment='left', verticalalignment='top',
         fontsize = 16, bbox=dict(facecolor='white', alpha=0.5), transform=axes.transAxes)
axes.set_xlim([10**28.7, 10**30.8])
axes.set_ylim([1e29, 10**35.5])
plt.minorticks_on()
plt.savefig(figurepath + 'localL-L.png')

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


# In[19]: Test density/luminosity functions
'''
for b in sdssfirst.bands:
    Z_cdf = np.arange(0.1, 5, 0.1)
    i = b.limited_indeces
    Z = sdssfirst.Z[i]
    L, Lmin = quas.localize(sdssfirst.Z[i], b.L[i], b.Lmin[i], b.k_g)
    
    plt.figure()
    plt.plot(sdssfirst.Z[i], np.log10(L), '.', markersize = 2)
    plt.plot(sdssfirst.Z[i], np.log10(Lmin))
    plt.title(r"$g(z) = \frac{(1 + z)^k}{1 + (\frac{1 + z}{Z_{cr}})^k}$, $Z_{cr} = 3.7$, $k_r = 4.82$")
    plt.xlabel(r"$z$")
    plt.ylabel(r"$\log(L_r')$")
    
    CDF = [quas.cdf(z, Z, L, Lmin) for z in Z_cdf]
    
    #log-log fit
    #    incr = [j for j in range(1, len(Z) - 1) if Z[j] != Z[j + 1]]
    Zfit = Z_cdf
    CDFfit = CDF
    fit = interp.UnivariateSpline(np.log10(1 + Zfit), np.log10(CDFfit), s = .2)
    dfit = fit.derivative()
    
    logZcurve = np.arange(0, 0.8, 0.02)
    Zcurve = 10**logZcurve - 1
    
    plt.figure()
    plt.plot(Z_cdf, CDF, '.')
    plt.plot(Zcurve, 10**fit(logZcurve))
    plt.title(r"Cumulative Density Function  for $L_{" + b.name + '}$')
    plt.minorticks_on()
    
    def devolution(z):
        Z = 1 + z
        logZ = np.log10(Z)
        sigma = 10**fit(logZ)
        return dfit(logZ) / quas.dV_dz(z) * sigma / Z
    
    rho = [devolution(z) for z in Z_cdf]
    plt.figure()
    plt.plot(Z_cdf, rho)
    plt.title(r"Density Evolution for $L_{" + b.name + '}$')
    plt.minorticks_on()

    LDF = [quas.ldf(Z_cdf[z], L, b.k_g, rho[z]) for z in range(len(Z_cdf))]
    plt.figure()
    plt.plot(Z_cdf, LDF)
    plt.title(r"Luminosity Density Function for $L_{" + b.name + '}$')
    plt.minorticks_on()
    
    L_lf = np.arange(min(np.log10(L)), max(np.log10(L)), 0.5)
    L_lf = 10**L_lf
    CLF = [quas.clf(l, Z, L, Lmin) for l in L_lf]
    plt.figure()
    plt.plot(np.log10(L_lf), np.log10(CLF), '.')
    plt.title(r"Cumulative Luminosity Function for $L_{" + b.name + '}$')
    plt.minorticks_on()
    
    LF = [quas.lf(l, L_lf, CLF) for l in L_lf]
    plt.figure()
    plt.plot(np.log10(L_lf), np.log10(LF), '.')
    plt.title(r"Luminosity Function for $L_{" + b.name + '}$')
    plt.minorticks_on()
'''
