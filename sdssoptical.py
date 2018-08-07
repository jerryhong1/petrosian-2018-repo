 # coding: utf-8

# In[1]: Pacakges
from __future__ import division
import numpy as np
#import scipy as sp
#from scipy.integrate import quad
import scipy.stats as stat
import math
import matplotlib.pyplot as plt
#import pandas as pd
from astropy.io import fits
from matplotlib import rc
import quasars as quas
from quasars import Band
from quasars import QuasarData
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# In[2]: upload .fit file
dataset = fits.open('dr7qso.fit')
dataset.info()

# In[3]: data labels for fits file:
#dataset[0].header

# In[4]: Access the data, get flux + redshift data

raw_data = dataset[1] 
#print quas.data[123] aze
Z = raw_data.data['z'] #redshifts
M = raw_data.data['IMAG'] #"photometric" magnitudes in this band
#print 'magnitudes: ', M

# In[6]: luminosity distance array — probably not necessary
dl = np.array([quas.d_lum(z) for z in Z])
#dl

# In[7]: convert magnitudes to fluxes (in cgs). see Abazajian et al. (2009) §3.3
# (multiply Jy by e-23 to get cgs: ergs/s/cm^2/Hz)

def magtoflux(m):
    f0 = 3.631e-20 
    return f0*10**(-m/2.5)

F = [magtoflux(m) for m in M]
F = np.array(F)

#convert to 2500 Angstrom band
alpha_opt = -0.5    
lambda_i = 7625. 
lambda_opt = 2500.
i_to_opt = lambda i: i *(lambda_opt / lambda_i)**(-alpha_opt)

F_opt = i_to_opt(F)

# In[8]: define k_opt

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
def k_opt(z): #with sign convention: m_intrinsic = m_obs - K, and L = 4 pi d^2 f/k(z)
    k_avg = -2.5 * (1 + alpha_opt) * math.log10(1 + z) 
    if (z > max(Z_k)):
        k = k_avg  #assume no emission line effect
    else:  
        k = np.interp(z, Z_k, K) - 2.5 * (1 + alpha_opt) * math.log10(1 + 2)
    return 10**(-k / 2.5)

fmin_i = 0.08317e-26 #see Singal et al. (2013)
fmin_opt = i_to_opt(fmin_i)

# In[11]: truncate data to obtain f_min = 0.083 mJy, aka m_max = 19.1
sdss_index = []
sdss_trunc_index = []
m_max = 19.1 
m_min = 15

for i in range(len(M)):
    if(M[i] > m_min):
        if(M[i] < m_max):
            sdss_index.append(i)
        else:
            sdss_trunc_index.append(i)

# In[9]: create objects
opt_band = Band('o', fmin_opt, F_opt[sdss_index], k_opt)
sdss = QuasarData(Z[sdss_index], [opt_band])

opt_band_trunc = Band('o', fmin_opt, F_opt[sdss_trunc_index], k_opt)
sdss_trunc = QuasarData(Z[sdss_trunc_index], [opt_band_trunc])

sdss.sort()
sdss_trunc.sort()

# In[12]: Test code for tau

i = 38123
j = [m for m in range(0, sdss.size()) if (opt_band.L[m] > opt_band.Lmin[i] and sdss.Z[m] < sdss.Z[i] and opt_band.Lmin[m] < opt_band.Lmin[i])]
j.append(i)
L_ass = opt_band.L[j]
Z_ass = sdss.Z[j]

L_rank = stat.rankdata(L_ass, 'max') #ranks of all the local luminosities
rank = L_rank[-1]  #associated set does not include data point itself, so -1 to avoid overcounting.
exp = 0.5 * (1 + len(L_ass))
resid = rank - exp
var = 1/12.0 * (len(L_ass)**2 - 1)
print '\ni = 0 tau test:'
print rank, exp, resid, var
print resid/math.sqrt(var)
print '\n'

# In[13]: make sanity-checking plot (log(L) vs. z)

#L vs z
plt.figure(1, figsize=(10,8))
plt.semilogy(sdss.Z, opt_band.L,'.', markersize=1, label="my data", color='black')
plt.plot(sdss_trunc.Z, opt_band_trunc.L,'.', markersize=1, label="truncated", color='#40E0D0')
#plt.plot(sdss.Z, np.log10(opt_band.Lmin),'-', linewidth = 2, label="min values", color='red')

#Lmin, Lmax vs z
z_graph = np.arange(0.01,5.5,0.05) 
Lmin_graph = [opt_band.min_luminosity(z) for z in z_graph]
Lmax_graph = [opt_band.luminosity(z, i_to_opt(magtoflux(15))) for z in z_graph]
plt.semilogy(z_graph, Lmin_graph,'-', markersize=3, label="max values", color='red')
plt.semilogy(z_graph, Lmax_graph,'-', markersize=3, label="max values", color='red')

#associated set (see above cell)
#plt.plot(Z_ass, np.log10(L_ass),'.', markersize=1, label="associated set for source", color='red')
#plt.plot(Z_ass[-1], np.log10(L_ass[-1]), '.', markersize=12, label="source", color = 'green')

#labeling
plt.xlabel("$z$", fontsize = 18)
plt.ylabel("$L_{opt}$ (erg s$^{-1}$ Hz$^{-1}$)", fontsize = 18)
plt.title("$L_{opt}$ at 2500 A vs.\ $z$ for SDSS DR7 Quasar Set", fontsize = 16)
#plt.legend(loc = "upper right")
axes = plt.gca()
axes.set_xlim([0,6])
axes.set_ylim([1e29,1e33])
plt.minorticks_on()
plt.savefig("../figures/SDSSlogLz.png")
plt.show()

# In[13.5] Plot demonstrating associated set
#L vs z
plt.figure(4, figsize=(10,8))
plt.semilogy(sdss.Z, opt_band.L,'.', markersize=1, label="my data", color='black')

#Lmin vs z
plt.semilogy(z_graph, Lmin_graph,'-', markersize=3, label="max values", color='red')

#associated set (see above cell)
lmin = opt_band.Lmin[i]
plt.semilogy(Z_ass, L_ass,'.', markersize=1, label="associated set for source", color='#900000')
plt.semilogy(Z_ass[-1], L_ass[-1], '.', markersize=12, label="source", color = 'red')
plt.semilogy([0,Z_ass[-1]], [lmin, lmin], linewidth = 1, color = 'red')
plt.semilogy([Z_ass[-1],Z_ass[-1]], [lmin ,1e40], linewidth = 1, color = 'red')

#labeling
plt.xlabel("$z$", fontsize = 18)
plt.ylabel("$L_{opt}$ (erg s$^{-1}$ Hz$^{-1}$)", fontsize = 18)
#plt.title("$\log(L)$ at 2500 A vs.\ $z$ for SDSS DR7 Quasar Set", fontsize = 16)
#plt.legend(loc = "upper right")
axes = plt.gca()
axes.text(0.5*(Z_ass[-1]), 10**(0.5*(np.log10(lmin) + 33)), 
        r'\begin{center} $\textbf{Associated}$ \\ $\textbf{Set}$ \end{center}',
        horizontalalignment='center',
        verticalalignment='center', 
        color = 'red', fontsize = 25)
axes.text(Z_ass[-1] + 0.07, opt_band.L[i], r'$\textbf{Source}$',
        horizontalalignment='left',
        verticalalignment='center', 
        color = 'white', fontsize = 20)

axes.set_xlim([0,5])
axes.set_ylim(np.array([1e29,1e33]))
plt.minorticks_on()
plt.savefig("../figures/AssociatedSet.png")
plt.show()

 # In[15]: Local Luminosity creation, given k:

#g(z) such that L' = L/g(z). Taken from Singal et al. (2013) and considers high-redshift objects.
def g(z, k):
    z_cr = 3.7
    return (1 + z)**k /(1 + ((1 + z) / z_cr)**k)

#test with k = 3:
k = 3.3
L_local, Lmin_local = quas.localize(sdss.Z, opt_band.L, opt_band.Lmin, k)
print 'L_local:', L_local, '\n'

#graph of L' vs. z
plt.figure(2)
plt.plot(sdss.Z, np.log10(L_local),'.', markersize=1, label="my data", color='black')
plt.title("log(L') vs.\ z for k = " + str(k))

plt.plot(sdss.Z, np.log10(Lmin_local), markersize=3, label="my data", color='red')
axes = plt.gca()
axes.set_xlim([0,6])
axes.set_ylim([28.5,32])

#tau(data[:,0], L_local, Lmin_local)

# In[17]: Iterate through k's to graph tau vs k (sample):

K = np.arange(3.1, 3.7, 0.1)
#K = np.array([3.3])
n = 1
Tau = quas.tauvsk(sdss.Z[::n], opt_band.L[::n], opt_band.Lmin[::n], K)

# In[18]: Tau vs k plot:
sigma = [-1, 1]
K = np.load('results/sdss-k.npy')
Tau = np.load('results/sdss-tau.npy')
k = np.interp(0, -Tau, K)
#tau vs k_opt
plt.figure(3)
plt.title(r"$\tau$ vs $k_{opt}$ (SDSS Optical Data)", fontsize = 14)
plt.plot(K, Tau)
axes = plt.gca()
axes.set_xlim([min(K), max(K)])
axes.set_ylim([-3,3])
plt.xlabel(r"$k_{opt}$", fontsize = 14)
plt.ylabel(r"$\tau$", fontsize = 14)

plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)

plt.plot([k, k], [-4, 0], color = 'red', linewidth = 1)
plt.text(k + 0.01, 0.1, r"$k_{opt} = $ " + str(round(k, 2)), color = 'red', fontsize = 14)
k_opterr = [np.interp(s, -Tau, K) for s in sigma]
plt.plot([k_opterr[0], k_opterr[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_opterr[1], k_opterr[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)

plt.savefig('../figures/tauvsk-sdss.eps')
plt.show()