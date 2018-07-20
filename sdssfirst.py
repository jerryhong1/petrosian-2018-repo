#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import math
import quasars as quas
import matplotlib.pyplot as plt
import scipy.stats as stat
#from mpl_toolkits.mplot3d import axes3d
from matplotlib import rc
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# In[1]: import file and put into array
file = open("SDSSFIRSTComplete.dat","r") 
lines = file.readlines()
file.close() 
data = np.empty((len(lines) - 5, 7))
for i in range(5, len(lines)):
    linearr = np.array(lines[i].strip().split())
    for j in range(0, len(linearr)):
        data[i - 5, j] = float(linearr[j])
    
# z = 0, d_lum = 1, k = 2, f_2500A = 3, f_1.4GHz = 4, L_o = 5, L_r = 6
    

# In[2]: Determine minimum values
fminOpt = 0.08318e-26 #0.083 mJy
alphaOpt = -0.5
def LMinOpt(z):
    dl = quas.d_lum(z)
    k = quas.k_opt(z)
    L = 4 * math.pi * dl**2 * fminOpt / k
    lambda_i = 7625.e-8 
    lambda_f = 2500.e-8
    L = L * (lambda_f / lambda_i)**(-alphaOpt)
    return L
    
alphaRad = -0.6
fminRad = 1.0e-26 #1 mJy
def LMinRad(z):
    dl = quas.d_lum(z)
    k = (1 + z)**(1 + alphaRad)
    L = 4 * math.pi * dl**2 * fminRad / k
    return L

# In[4]: Localize test
Lmin_opt = np.array([LMinOpt(z) for z in data[:,0]])
Lmin_rad = np.array([LMinRad(z) for z in data[:,0]])
L_local, Lmin_local = quas.localize(data[:,0], data[:,5], Lmin_opt, 3.5)

#plt.figure()
#plt.plot(data[:,0], np.log10(L_local),'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10(Lmin_local))
#print quas.tau(data[:,0], L_local, Lmin_local)

# In[5]: SKIP combined tau analysis
'''
def tauVsKs(Z, L_opt, Lmin_opt, L_rad, Lmin_rad, K_opt, K_rad):
    Opt = []
    Rad = []
    Tau = []
    for ko in K_opt:
        for kr in K_rad:
            L_opt_l, Lmin_opt_l = quas.localize(Z, L_opt, Lmin_opt, ko)
            L_rad_l, Lmin_rad_l = quas.localize(Z, L_rad, Lmin_rad, kr)
            print('\nko = ' +  str(ko) + '  kr = ' + str(kr))
            t = tauCombined(Z, L_opt_l, Lmin_opt_l, L_rad_l, Lmin_rad_l)
            print ('\ntau = ' + str(t))
            Opt.append(ko)
            Rad.append(kr)
            Tau.append(t)
    return Opt, Rad, Tau

#modified defintion of associated set to incorporate both bands:
def tauCombined(Z, L_opt, Lmin_opt, L_rad, Lmin_rad):
    #as defined in singal, petrosian papers in 2010s. tau = (∑resid)/(sqrt(∑variance))
    resido = 0
    varo = 0
    
    residr = 0
    varr = 0
    
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L_opt[m] > Lmin_opt[i] and L_rad[m] > Lmin_rad[i] and Z[m] < Z[i])] #see petrosian
        if (len(j) == 1 or len(j) == 2): continue
        j.append(i)
        
        #determine rank: optical
        L_ass_opt =  L_opt[j]
        L_rank = stat.rankdata(L_ass_opt, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (len(L_ass_opt))
        
        resido = resido + (rank - exp)
        varo = varo + (1/12.0 * (-1 + (len(L_ass_opt)-1)**2))

        #determine rank: optical
        L_ass_rad =  L_rad[j]
        L_rank = stat.rankdata(L_ass_rad, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (len(L_ass_rad))
        
        residr = residr + (rank - exp)
        varr = varr + (1/12.0 * (-1 + (len(L_ass_rad)-1)**2))
        
        #troubleshooting
        if(i % 500 == 0): print i, resido, varo, residr, varr
        
    t = residr**2 / varr +  resido**2 / varo
    return math.sqrt(t)

o, r, t = tauVsKs(data[:,0], data[:,5], Lmin_opt, data[:,6], Lmin_rad, np.arange(3.0, 3.4, 0.1), np.arange(5.3, 5.7, 0.1))
print o, r, t


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot a basic wireframe.
ax.plot_wireframe(o, r, t, rstride=10, cstride=10)
plt.show()
'''

# In[6]: Determine if "optically limited" or "radio limited"
Z = np.hstack((np.arange(0.001, 4.99, 0.001), np.arange(5, 500, 0.5)))
o = lambda logl: np.interp(logl, [np.log10(LMinOpt(z)) for z in Z], Z, 0, float("inf"))
ZMaxOpt = o(np.log10(data[:,5]))

r = lambda logl: np.interp(logl, [np.log10(LMinRad(z)) for z in Z], Z, 0, float("inf"))
ZMaxRad  = r(np.log10(data[:,6]))

#data points that are optically limited, data points that are radio limited:
OptLim = np.where(ZMaxOpt < ZMaxRad)[0]
RadLim = np.where(ZMaxRad < ZMaxOpt)[0]

plt.figure()
plt.plot(data[OptLim,0], np.log10(data[OptLim,5]),'.', markersize=1, label="my data", color='black')
plt.plot(data[:,0], np.log10([LMinOpt(z) for z in data[:,0]]))

plt.figure()
plt.plot(data[RadLim,0], np.log10(data[RadLim,6]),'.', markersize=1, label="my data", color='black')
plt.plot(data[:,0], np.log10([LMinRad(z) for z in data[:,0]]))


K_opt = np.arange(2.5, 4., 0.1)
K_rad = np.arange(4.5, 5.5, 0.1)
Tau_opt = quas.tauvsk(data[OptLim,0], data[OptLim,5], Lmin_opt[OptLim], K_opt)
Tau_rad = quas.tauvsk(data[RadLim,0], data[RadLim,6], Lmin_rad[RadLim], K_rad)

# In[7]: Corr-reduced luminosity and correlations

#alpha = 1 #test
L0 = 1e30 #test
k_opt = np.interp(0, -Tau_opt, K_opt)
k_rad = np.interp(0, -Tau_rad, K_rad)

L_optl, foo = quas.localize(data[:,0], data[:,5], Lmin_opt, k_opt)
L_radl, foo = quas.localize(data[:,0], data[:,6], Lmin_rad, k_rad)

LumCr = lambda L_opt, L_rad, alpha: L_rad * (L0/L_opt)**alpha
#LCr = [LumCr(L_optl[j], L_radl[j], alpha) for j in range(0, len(data[:,]))]

def rvsalpha(L_optl, L_radl, Alpha):
    R = []
    for a in Alpha:
        LCr = [LumCr(L_optl[i], L_radl[i], a) for i in range(0, len(data[:,]))]        
        r = stat.linregress(np.log10(L_optl), np.log10(LCr)).rvalue
        R.append(r)
    return R

Alpha = np.arange(0,1,0.005)
R = rvsalpha(L_optl, L_radl, Alpha)


# In[3]: Plots of log(L) vs. z
plt.figure()
plt.plot(np.log10(data[:,5]), np.log10(data[:,6]), '.', markersize=1, label="my data", color='black')
plt.xlabel(r"$\log(L_{opt})$")
plt.ylabel(r"$\log(L_{rad})$")
plt.title("Optical vs. Radio")
plt.minorticks_on()

#logl_opt vs z
plt.figure(figsize=(8,6))
plt.plot(data[:,0], np.log10(data[:,5]),'.', markersize=1, label="my data", color='black')
#plt.plot(data[:,0], np.log10([LMinOpt(z) for z in data[:,0]]), linewidth = 1)
plt.xlabel(r"$z$", fontsize=14)
plt.ylabel(r"$\log(L_{opt})$", fontsize=14)
plt.title(r"$\log(L)$ at 2500 \r{A} vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([29, 33])
plt.minorticks_on()
plt.savefig('../figures/OptSDSSFIRSTlogLz.png')

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
plt.ylabel(r"$\log(L_{rad})$", fontsize=14)
plt.title(r"$\log(L)$ at 1.4GHz vs.\ $z$ for SDSS X FIRST Set", fontsize=16)
plt.minorticks_on()
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([29, 36.5])
plt.savefig('../figures/RadSDSSFIRSTlogLz.png')
plt.show()

'''
#Test Zmax
plt.plot(data[i,0], np.log10(data[i,6]), '.', markersize=12, color = 'red')
plt.plot(ZMaxRad[i], np.log10(data[i,6]), '.', markersize=12, color = 'red')
'''
# In[] Zmax Definition

plt.figure(figsize=(8,6))
#logl_rad vs z
plt.plot(data[:,0], np.log10(data[:,6]),'.', markersize=1, label="my data", color='#A0A0A0')
plt.plot(data[:,0], np.log10([LMinRad(z) for z in data[:,0]]), linewidth = 1, color = 'red')
plt.xlabel("z", fontsize=14)
plt.ylabel(r"$\log(L_{rad})$", fontsize=14)
plt.minorticks_on()
axes = plt.gca()
axes.set_xlim([0, 5])
axes.set_ylim([29, 36.5])

#Zmax
i = 2002
plt.plot(data[i,0], np.log10(data[i,6]), '.', markersize=12, color = 'black')
plt.plot([data[i,0],ZMaxRad[i]], [np.log10(data[i,6]), np.log10(data[i,6])], '--', color = 'black', linewidth = 2)
plt.plot([ZMaxRad[i],ZMaxRad[i]], [0, np.log10(data[i,6])], '--', color = 'black', linewidth = 2)
axes.text(data[i,0], np.log10(data[i,6]) + 0.3, r'\textbf{Source $\mathbf{i}$}',
        horizontalalignment='center',
        verticalalignment='center', 
        color = 'black', fontsize = 14,
        bbox=dict(facecolor='white', alpha=0.7, ec = 'none'))
axes.text(ZMaxRad[i] + 0.1, 30, r'$z_{max, i}^{r}$',
        horizontalalignment='left',
        verticalalignment='top', 
        color = 'black', fontsize = 18)
plt.savefig('../figures/zmax.png')
plt.show()

# In[] Plots of tau vs k
sigma = [-1, 1]
#tau vs k_opt
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
plt.text(k_opt + 0.01, 0.1, r"$k_{opt} = $ " + str(round(k_opt, 2)), color = 'red', fontsize = 14)
k_opterr = [np.interp(i, -Tau_opt, K_opt) for i in sigma]
plt.plot([k_opterr[0], k_opterr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_opterr[1], k_opterr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)

plt.savefig('../figures/tauvsk-first-opt.eps')
plt.show()

#tau vs k_rad
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
plt.text(k_rad + 0.01, 0.1, r"$k_{rad} = $ " + str(round(k_rad, 2)), color = 'red', fontsize = 14)
k_raderr = [np.interp(i, -Tau_rad, K_rad) for i in sigma]
plt.plot([k_raderr[0], k_raderr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_raderr[1], k_raderr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)
         
plt.savefig('../figures/tauvsk-first-rad.eps')
plt.show()

# In[8]: More plots
#Lopt' vs Lrad' (local)
plt.figure(figsize=(8,6))
plt.plot(np.log10(L_optl), np.log10(L_radl), '.', markersize=1)
plt.xlabel(r"$\log(L_{opt}')$", fontsize = 16)
plt.ylabel(r"$\log(L_{rad}')$", fontsize = 16)
#plt.title("Local luminosity scatter plot")
plt.text(29.75, 35, r"$g(z) = \frac{(1 + z)^k}{1 + (\frac{1 + z}{Z_{cr}})^k}$", color = 'black', fontsize = 18)
plt.text(30.5, 34.5, r"$k_{rad} = $ " + str(round(k_rad, 2)) + "\n" + r"$k_{opt} = $ " + str(round(k_opt, 2)) + "\n" + r"$Z_{cr} = $ " + str(3.7), color = 'black', fontsize = 16, bbox=dict(facecolor='white', alpha=0.5))
axes = plt.gca()
axes.set_xlim([28.5, 31])
axes.set_ylim([29, 36])
plt.savefig('../figures/localL-L.png')
plt.minorticks_on()

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
plt.savefig('../figures/r-alpha.eps')