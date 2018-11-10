'''
Author: Jerry Hong. 
Aug 27, 2018

Professor Singal suggested that I try different associated sets and see what I 
get for tau vs k. I will try that on a gamma-limited (when comparing to optical)
dataset. 

1. Gamma-ray: 297 FSRQs, Gamma = 2.45 assumed.
2. Optical: 4211 quasars from SDSSxFIRST  set
3. Radio: 1234 quasars from SDSSxFIRST set

universal flux limit, Z_cr = 3.7 assumed. All data is sorted by redshift, which is important.

There is a variable called type:
    1: associated set is all j s.t. z_j < z_i, L_j > Lmin_i. This is what I have been using.
    2: associated set is all j s.t. z_j < z_i, L_j > Lmin_i, zmax_j < zmax_i.
    3: "flattens out" the Lmin curve and proceeds as type 1.
'''

# In[1]: import data; will not consolidate in an QuasarData object, since it's just
# one band and one
import quasars as quas
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as stat
from scipy.signal import argrelextrema, find_peaks
plt.style.use('mystyle.mplstyle')
log = np.log10
# 'sdssfirst' or 'fsrq'
data = 'fsrq'

if data == 'fsrq':
    # just gamma-limited fsrq's
    L = np.load('results/fsrq/L_gamma.npy')
    Z = np.load('results/fsrq/Z_gamma.npy')
    Lmin = np.load('results/fsrq/Lmin_gamma.npy')
    
    fmin = 3e-12
    Gamma = 2.44
    alpha = 1 - Gamma
    kcorrect = lambda z: quas.kcorrect(z, alpha)
    
    ind = 240
    k = 6
    K_arr = np.arange(4.5, 6.5, 0.1)

elif data == 'sdssfirstopt':
    # just gamma-limited fsrq's
    L = np.load('results/sdssfirst/L_opt.npy')
    Z = np.load('results/sdssfirst/Z_opt.npy')
    Lmin = np.load('results/sdssfirst/Lmin_opt.npy')
    
    fmin = quas.bandtoband(0.08317e-26, quas.lambda_i, quas.lambda_opt, quas.alpha_opt)
    kcorrect = lambda z: quas.kcorrect_SDSS(z)
    ind = 2604
    k = 4
    K_arr = np.arange(3, 4, 0.1)

elif data == 'sdssfirstrad':
    # just gamma-limited fsrq's
    L = np.load('results/sdssfirst/L_rad.npy')
    Z = np.load('results/sdssfirst/Z_rad.npy')
    Lmin = np.load('results/sdssfirst/Lmin_rad.npy')
    
    fmin = 1.0e-26 #1 mJy
    kcorrect = lambda z: quas.kcorrect_rad(z)

    ind = 1054
    k = 5
    K_arr = np.arange(4, 6, 0.1)

def luminosity(z, f): 
    kc = kcorrect(z)
    return 4 * math.pi * (quas.d_lum(z)**2) * f / kc

# In[3]: tau, flatten functions
def tau(Z, L, Lmin, Zmax, type): 
    resid = 0
    var = 0
    
    if type == 3:
        Lmin = flatten(Z, Lmin)
        t = [i for i in range(len(Z)) if L[i] >= Lmin[i]]
        Z = Z[t]
        L = L[t]
        Lmin = Lmin[t]
        Zmax = Zmax[t]
        print len(Z)
        
    for i in range(len(Z)):
        #create associated sets
        if type == 1 or type == 3:
            j = [m for m in range(len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])]
        elif type == 2:
            j = [m for m in range(len(Z)) if (Zmax[m] > Z[i] and Z[m] < Z[i])]
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


# assumes Lmin hits max and dips back down, and then returns up. assume data is sorted.
def flatten(Z, Lmin):
    i_1 = find_peaks(Lminl)[0]
    if i_1.size == 0:
        # if there is no max, monotonicity assumed.
        return Lmin
    else: 
        i_1 = i_1[0]
    Lmin_max = Lmin[i_1]
    Z_1 = Z[i_1]
    i_2 = [i for i in range(len(Z)) if (Z[i] > Z_1 and Lmin[i] > Lmin_max)]
    if i_2:
        Z_2 = Z[i_2[0]]
    else: Z_2 = float('inf')    
    
    return np.array([Lmin_max if Z[i] < Z_2 and Z[i] > Z_1 else Lmin[i] for i in range(len(Z))])

# assume data is sorted.
def Zmax(L, luminosity, fmin, kcorrect, k, Z): 
    Z_zmax = np.hstack((np.arange(0.01, 4.99, 0.01), np.arange(5, 500, 0.5)))
    log_Lmin = np.array([log(luminosity(z, fmin) / quas.g(z, k)) for z in Z_zmax])
#    plt.semilogy(Z_zmax[5:500], log_Lmin[5:500], '.')
    i_max = find_peaks(log_Lmin)[0] #see if there is a local max
    if i_max.size > 0:
        i_max = i_max[0]
    else:
        i_max = -1

    def interp(i):
        l = L[i]
        if l== 0.0: return 0.0 
        if i_max >= 0:
            if Z[i] < Z_zmax[i_max] and log(L[i]) < log_Lmin[i_max]:
                return np.interp(log(l), log_Lmin[:i_max], Z_zmax[:i_max], 0, float("inf"))
        return np.interp(log(l), log_Lmin, Z_zmax, 0, float("inf"))
    
    return np.array([interp(i) for i in range(len(L))])


def tauvsk(Z, L, Lmin, Zmax, K, type):
    s = []
    for k in K:
        Ll, Lminl = quas.localize(Z, L, Lmin, k)
        Zmaxl = Zmax(Ll, luminosity, fmin, kcorrect, k, Z)
        print('\nk = ' +  str(k))
        t = tau(Z, Ll, Lminl, Zmaxl, type)
        print ('\ntau = ' + str(t))
        s.append(t)
    
    s = np.array(s)
    print 'tau = 0 at k = ' + str(np.interp(0, -s, K))
    return np.array(s)

# In[2]: visualize associated sets (k = 5.5)
Ll, Lminl = quas.localize(Z, L, Lmin, k)

Zmaxl = Zmax(Ll, luminosity, fmin, kcorrect, k, Z)

#type = 1:
j = [m for m in range(len(Z)) if (Ll[m] > Lminl[ind] and Z[m] < Z[ind])]
plt.figure()
plt.semilogy(Z, Ll, 'k.', Z, Lminl, 'b')
plt.semilogy(Z[ind], Ll[ind], 'x', markersize = 15)
plt.semilogy(Z[j], Ll[j], '.', color = 'red')


#type = 2:
j = [m for m in range(len(Z)) if (Zmaxl[m] > Z[ind] and Z[m] < Z[ind])]
plt.figure()
plt.semilogy(Z, Ll, 'k.', Z, Lminl, 'b')
plt.semilogy(Z[ind], Ll[ind], 'x', markersize = 15)
plt.semilogy(Z[j], Ll[j], '.', color = 'red')


#type = 3: 
Lminl_flat = flatten(Z, Lminl)
j = [m for m in range(len(Z)) if (Ll[m] > Lminl_flat[ind] and Z[m] < Z[ind])]
t = [i for i in range(len(Z)) if Ll[i] >= Lminl_flat[i]]
tnot = [i for i in range(len(Z)) if Ll[i] < Lminl_flat[i]]
plt.figure()
plt.semilogy(Z[tnot], Ll[tnot], '.', alpha = 0.5, color = '#909090') 
plt.semilogy(Z, Lminl_flat)
plt.semilogy(Z[t], Ll[t], '.', color = 'black')
plt.semilogy(Z[ind], Ll[ind], 'x', markersize = 15)
plt.semilogy(Z[j], Ll[j], '.', color = 'red')

# In[]: Calculate tau:
Tau_arr1 = tauvsk(Z, L, Lmin, Zmax, K_arr, 1)
Tau_arr2 = tauvsk(Z, L, Lmin, Zmax, K_arr, 2)
Tau_arr3 = tauvsk(Z, L, Lmin, Zmax, K_arr, 3)

# In[]: Plot tau vs k:

k = np.interp(0, -Tau_arr1, K_arr)
K = K_arr
sigma = [-1, 1]
plt.figure()
for Tau in [Tau_arr1, Tau_arr2, Tau_arr3]:
    plt.plot(K, Tau)
axes = plt.gca()
axes.set_xlim([min(K), max(K)])
axes.set_ylim([-3,3])
plt.xlabel(r"$k_{\gamma}$")
plt.ylabel(r"$\tau$")

plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)

'''
plt.plot([k, k], [-4, 0], color = 'red', linewidth = 1)
k_err = [np.interp(i, -Tau, K) for i in sigma]
plt.text(k + 0.01, 0.1, r"$k_{\gamma} = " + str(round(k, 2)) 
    + '^{ + ' + str(round(k_err[1] - k, 2)) + '}_{' + str(round(k_err[0] - k, 2)) + '}$' ,
    color = 'red', fontsize = 24)
plt.plot([k_err[0], k_err[0]], [-4, 1], '--', color = 'red', linewidth = 0.5)
plt.plot([k_err[1], k_err[1]], [-4, -1], '--', color = 'red', linewidth = 0.5)
'''