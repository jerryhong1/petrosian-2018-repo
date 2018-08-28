#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import quasars as quas
from quasars import Band, QuasarData
import scipy.interpolate as interp
import matplotlib.pyplot as plt
plt.style.use('mystyle.mplstyle')

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

# In[3]: Create Objects
o_band = Band('opt', fmin_opt, data[:,3], quas.kcorrect_SDSS)
o_band.set_luminosity(data[:,5])
r_band = Band('rad', fmin_rad, data[:,4], quas.kcorrect_rad)
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
figurepath = '../figures/sdssxfirst/'
plt.figure()
plt.loglog(data[:,5], data[:,6], '.', markersize=1, label="my data", color='black')
plt.xlabel(r"$\log(L_{opt})$")
plt.ylabel(r"$\log(L_{rad})$")
plt.title("Optical vs. Radio")

#logl_opt vs z
plt.figure()
plt.semilogy(data[:,0], data[:,5],'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10([LMinOpt(z) for z in data[:,0]]), linewidth = 1)
plt.xlabel(r"$z$", fontsize=14)
plt.ylabel(r"$L_{opt}$", fontsize=14)
plt.title(r"Luminosity at 2500 \r{A} vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([1e29, 1e33])
plt.savefig(figurepath + 'OptSDSSFIRSTlogLz.png')

#logl_rad vs z
plt.figure()
plt.plot(data[:,0], np.log10(data[:,6]),'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10([LMinRad(z) for z in data[:,0]]), linewidth = 1)
plt.xlabel(r"$z$", fontsize=14)
plt.ylabel(r"$L_{rad}$", fontsize=14)
plt.title(r"Luminosity at 1.4GHz vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([29, 36.5])
plt.savefig(figurepath + 'RadSDSSFIRSTlogLz.png')
plt.show()

# In[] PLOT: Zmax + Associated Set Definition
plt.figure(figsize=(8,6))
#logl_rad vs z
plt.semilogy(data[:,0], r_band.L,'.', markersize=0.5, label="my data", color='black')
plt.semilogy(data[:,0], [r_band.min_luminosity(z) for z in data[:,0]], 
             linewidth = 1, color = 'red')
plt.xlabel("z", fontsize=14)
plt.ylabel(r"$L_{rad}$", fontsize=14)
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([1e29, 3.16e36])

i = 3060

# associated set
j = [m for m in range(sdssfirst.size()) if (r_band.L[m] > r_band.Lmin[i] 
    and sdssfirst.Z[m] < sdssfirst.Z[i] and r_band.Lmin[m] < r_band.Lmin[i])]
j.append(i)
lmin = r_band.Lmin[i]
Z_ass = sdssfirst.Z[j]
L_ass = r_band.L[j]
plt.semilogy(Z_ass, L_ass,'.', markersize=0.5, label="associated set for source", color='#FF9999')
plt.semilogy([0,Z_ass[-1]], [lmin, lmin], linewidth = 1, color = 'red')
plt.semilogy([Z_ass[-1],Z_ass[-1]], [lmin ,1e40], linewidth = 1, color = 'red')
axes.text(0.5*(Z_ass[-1]), 10**(0.5*(np.log10(lmin) + 37)), 
        r'\begin{center} $\textbf{Associated}$ \\ $\textbf{Set}$ \end{center}',
        horizontalalignment='center',
        verticalalignment='center', 
        color = '#FF2222', fontsize = 18)

# Zmax    
plt.semilogy(data[i,0], data[i,6], '.', markersize=12, color = 'red')
plt.semilogy([0,r_band.Zmax[i]], [r_band.L[i], r_band.L[i]],
             '--', color = 'red', linewidth = 1)
plt.semilogy([r_band.Zmax[i],r_band.Zmax[i]], [1, r_band.L[i]], '--', 
             color = 'red', linewidth = 1)
axes.text(data[i,0], data[i,6] * 10**0.3, r'\textbf{Source} $i$',
        horizontalalignment='center',
        verticalalignment='center', 
        color = 'black', fontsize = 14,
        bbox=dict(facecolor='white', alpha=0.8, ec = 'none'))
axes.text(r_band.Zmax[i] + 0.05, 3e29, r'$z_{max, i}$',
        horizontalalignment='left',
        verticalalignment='top', 
        color = 'black', fontsize = 18)
axes.text(0.05, r_band.Lmin[i], r'$L_{min, i}$',
        horizontalalignment='left',
        verticalalignment='bottom', 
        color = 'black', fontsize = 14,
        bbox=dict(facecolor='white', alpha=0.7, ec = 'none'))
axes.text(0.05, r_band.L[i], r'$L_{i}$',
        horizontalalignment='left',
        verticalalignment='bottom', 
        color = 'black', fontsize = 14,
        bbox=dict(facecolor='white', alpha=0.7, ec = 'none'))
plt.savefig('../figures/associated-set-zmax.png')
plt.show()

# In[] PLOT of tau vs k
sigma = [-1, 1]
#tau vs k_opt
for b in sdssfirst.bands:
    Tau = b.tau_array
    K = b.k_array
    k = b.k_g
    
    plt.figure(3)
    plt.title(r"$\tau$ vs $k_{" + b.name + "}$ (SDSS X FIRST set)")
    plt.plot(K, Tau)
    axes = plt.gca()
    axes.set_xlim([min(K), max(K)])
    axes.set_ylim([-3,3])
    plt.xlabel(r"$k_{" + b.name + "}$", fontsize = 14)
    plt.ylabel(r"$\tau$", fontsize = 14)
    
    plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
    plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
    plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)
    
    plt.plot([k, k], [-4, 0], color = 'red', linewidth = 1)
    k_err = [np.interp(i, -Tau, K) for i in sigma]
    plt.text(k + 0.01, 0.1, r"$k_{" + b.name + "} = " + str(round(k, 2)) 
        + '^{ + ' + str(round(k_err[1] - k, 2)) + '}_{' + str(round(k_err[0] - k, 2)) + '}$' ,
        color = 'red', fontsize = 24)
    plt.plot([k_err[0], k_err[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
    plt.plot([k_err[1], k_err[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)
    
    plt.savefig(figurepath + 'tauvsk-first-' + b.name + '.eps')
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
plt.text(.53, 0.96, r"$k_{rad} = $ " + str(round(r_band.k_g, 2)) + "\n" + r"$k_{opt} = $ " 
         + str(round(o_band.k_g, 2)) + "\n" + r"$Z_{cr} = $ " + str(3.7), color = 'black',
         horizontalalignment='left', verticalalignment='top',
         fontsize = 16, bbox=dict(facecolor='white', alpha=0.5), transform=axes.transAxes)
axes.set_xlim([10**28.7, 10**30.8])
axes.set_ylim([1e29, 10**35.5])
plt.savefig(figurepath + 'localL-L.png')

#alpha vs r
plt.figure()
plt.plot(Alpha, R, color = 'black', linewidth = 1)
#plt.plot(Alpha[::10], R[::10], marker = '.', markersize = 8)
alpha = np.interp(0, -np.array(R), Alpha)
plt.plot([alpha, alpha], [-4, 0], color = 'red', linewidth = 1)
plt.ylabel(r"Correlation Coefficient")
plt.xlabel(r"$\alpha$")
#plt.title(r"Correlation of $L_{cr}^{\prime\; rad} = L_{rad}'(L_0/L_{opt}')^\alpha$ vs. $L_{opt}'$")
plt.title(r"Optical-Radio Correlation", fontsize = 16) # for poster
plt.text(alpha + 0.02, 0.01, r"$\alpha = $ " + str(round(alpha, 2)), color = 'red', fontsize = 14)

axes = plt.gca()
plt.plot(np.arange(-1, 3), np.zeros(4), color = 'black', linewidth = 1)
axes.set_xlim([min(Alpha), 1])
axes.set_ylim([min(R), max(R)])
plt.savefig(figurepath + 'r-alpha.eps')


# In[19]: Test density/luminosity functions

b = o_band
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

# In[20]: CDF fit, plot
#log-log fit
#    incr = [j for j in range(1, len(Z) - 1) if Z[j] != Z[j + 1]]
Zfit = Z_cdf
CDFfit = CDF
w = np.hstack((np.zeros(0) + 3., np.zeros(5) + .5, np.zeros(5) + .5, np.arange(len(Zfit) - 10) + 2)) # first 10 at weight 0.5 for optical
fit = interp.UnivariateSpline(np.log10(1 + Zfit), np.log10(CDFfit), s = .3, w = w)
dfit = fit.derivative()

logZcurve = np.arange(0, 0.8, 0.02)
Zcurve = 10**logZcurve - 1

plt.figure()
plt.plot(Z_cdf, CDF, '.')
plt.plot(Zcurve, 10**fit(logZcurve))
plt.title(r"Cumulative Density Function  for $L_{" + b.name + '}$')

def devolution(z):
    Z = 1 + z
    logZ = np.log10(Z)
    sigma = 10**fit(logZ)
    return dfit(logZ) / quas.dV_dz(z) * sigma / Z

rho = [devolution(z) for z in Z_cdf]
rho = np.array(rho) * 3.086e+24**3 #convert from /cm^3 to /Mpc^3
plt.figure()
plt.plot(Z_cdf, rho)
plt.title(r"Density Evolution for $L_{" + b.name + '}$')

LDF = [quas.ldf(Z_cdf[z], L, b.k_g, rho[z]) for z in range(len(Z_cdf))]
plt.figure()
plt.plot(Z_cdf, LDF)
plt.title(r"Luminosity Density Function for $L_{" + b.name + '}$')

# In[] LF

L_lf = np.arange(min(np.log10(L)), max(np.log10(L)), 0.5)
L_lf = 10**L_lf
CLF = [quas.clf(l, Z, L, Lmin) for l in L_lf]
plt.figure()
plt.plot(np.log10(L_lf), np.log10(CLF), '.')
plt.title(r"Cumulative Luminosity Function for $L_{" + b.name + '}$')

LF = [quas.lf(l, L_lf, CLF) for l in L_lf]
plt.figure()
plt.plot(np.log10(L_lf), np.log10(LF), '.')
plt.title(r"Luminosity Function for $L_{" + b.name + '}$')