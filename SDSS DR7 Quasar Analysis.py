 
# coding: utf-8

# In[1]: Pacakges
import numpy as np
import scipy as sp
from scipy.integrate import quad
import scipy.stats as stat
import math
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits

# In[2]: upload .fit file
dataset = fits.open('dr7qso.fit')
dataset.info()

# In[3]: data labels for fits file:
#dataset[0].header

# In[4]: Access the data, get flux + redshift data

quas = dataset[1]  #short for quasars
#print quas.data[123] 
Z = quas.data['z'] #redshifts
M = quas.data['ITMAG'] #"photometric" magnitudes in this band
#print 'magnitudes: ', M

# In[5]: luminosity distance in cm. assume the standard ΛCDM cosmology: H0 = 71 km s−1 Mpc−1, ΩΛ = 0.7, and Ωm = 0.3
#work in cgs
def integrand(x, a, b):
    return 1 / math.sqrt(a + b*(1+x)**3)

def d_lum(z):
    Omega_l = 0.7
    Omega_m = 0.3
    cH0 = 1.30291e+28 #cm, ratio of c/H0. H0 = 71 km/s/Mpc
    integral = quad(integrand, 0, z, args = (Omega_l, Omega_m))
    return integral[0] * (1+z) * cH0


# In[6]: luminosity distance array — probably not necessary
dl = np.array([d_lum(z) for z in Z])
#dl

# In[7]: convert magnitudes to fluxes (in cgs). see Abazajian et al. (2009) §3.3
# (multiply Jy by e-23 to get cgs: ergs/s/cm^2/Hz)

def magtoflux(m):
    f0 = 3.631e-20 
    return f0*10**(-m/2.5)

F = [magtoflux(m) for m in M]
F = np.array(F)


# In[8]: calculate luminosity, minimum luminosity

#k-corrections for richard et al. (2006)
f = open('kcorr_quas.txt', 'r')
K = []
Z_k = []
for line in f:
    Z_k.append(line[0:4])
    K.append(line[5:11])

#print K, Z_k

def k(z): #with sign convention: m_intrinsic = m_obs - K, and L = 4 pi d^2 f/k(z)
    k = float(K[int(100 * z)])            
    return 10**(-k / 2.5)

def lum(z, f):
    alpha = -0.5 #for optical k-correction
#    kc = (1+z)**alpha
    kc = k(z)
    
    #convert to 2500 Angstrom band
    lambda_i = 7625. 
    lambda_f = 2500.
    f_adj = f*(lambda_f / lambda_i)**(-alpha)
    return 4*math.pi*(d_lum(z)**2)*f_adj/kc

fmin = 0.083e-26 #see Singal et al. (2013)
def lmin(l, f):
    return l*fmin/f

L = [lum(Z[i], F[i]) for i in range(0,len(F))]
L = np.array(L)
Lmin = [lmin(L[i], F[i]) for i in range(0, len(F))]
Lmin = np.array(Lmin)

#print(np.argmax(L))
#L, Lmin
# In[9]: Columns of data:
# 0 = z
# 
# 1 = M
# 
# 2 = F
# 
# 3 = L
# 
# 4 = Lmin

# In[10]: concatenate all data and calcs
data = np.stack((Z, M, F, L, Lmin))
data = np.transpose(data)
df = pd.DataFrame(data, columns= ['z', 'M', 'F', 'L', 'Lmin'])

# In[11]: truncate data to obtain f_min = 0.083 mJy, aka m_max = 19.1
data_temp = [];
data_trunc = [];
m_max = 19.1 
m_min = 15
for i in range(0, len(data[:,1])):
    if(data[i,1] > m_min):
        if(data[i,1] < m_max):
            data_temp.append(data[i,:])
        else:
            data_trunc.append(data[i,:])

data_trunc = np.array(data_trunc)        
data = np.array(data_temp)
print data[:,1], data.shape, np.amax(data[:,0])


# In[12]: Test code for tau; probably not necessary.

i = 8124
#k = 0
j = [m for m in range(0, len(data[:,])) if (data[m,3] > data[i,4] and data[m,0] < data[i,0])]
#j = [m for m in range(0, len(data[:,])) if (data[m,3] >= data[i,4])]
j.append(i)
#rankIndex = np.where(np.array(j) == i)
L_ass = data[j,3]
Z_ass = data[j,0]

#L_local = [L_ass[p] * (1 + Z_ass[p])**(-k) for p in range(0, len(L_ass))]
L_rank = stat.rankdata(L_ass, 'max') #ranks of all the local luminosities
rank = L_rank[-1] - 1 #associated set does not include data point itself, so -1 to avoid overcounting.
exp = 0.5 * (1 + len(L_ass))
resid = rank - exp
var = 1/12.0 * (len(L_ass)**2 - 1)
print '\ni = 0 tau test:'
print rank, exp, resid, var
print resid/math.sqrt(var)
print '\n'

# In[13]: make sanity-checking plots

#L vs z
plt.figure(1, figsize=(10,8))
plt.plot(data[:,0], np.log10(data[:,3]),'.', markersize=1, label="my data", color='black')
plt.plot(data_trunc[:,0], np.log10(data_trunc[:,3]),'.', markersize=1, label="truncated", color='#40E0D0')

#Lmin vs z
zmin = data[:,0][np.argsort(data[:,0])]
Lmin_sorted = data[:,4][np.argsort(data[:,0])]
plt.plot(zmin, np.log10(Lmin_sorted),'-', markersize=3, label="min values", color='red')

#Lmax vs z

Lmax = [lum(z, magtoflux(15)) for z in np.arange(0.01,5,0.05)]
plt.plot(np.arange(0.01,5,0.05), np.log10(Lmax),'-', markersize=3, label="max values", color='red')

#associated set (see above cell)
plt.plot(Z_ass, np.log10(L_ass),'.', markersize=1, label="associated set for source", color='red')
plt.plot(Z_ass[-1], np.log10(L_ass[-1]), '.', markersize=12, label="source", color = 'green')

#labeling
plt.xlabel("z")
plt.ylabel("log(L)")
plt.title("log(L) at 2500 A vs. z for SDSS DR7 Quasar Set")
#plt.legend(loc = "upper right")
axes = plt.gca()
axes.set_xlim([0,6])
axes.set_ylim([29,33])
plt.show()


# In[14]: kendall's tau statistic: takes in L'-z coordinates (+Lmin') to and calculates correlation, 
# constructing associated sets for each source along the way. Assume Z, L, Lmin are all same length.
def tau_(Z, L, Lmin): #as defined in singal, petrosian papers in 2010s. tau = (∑resid)/(sqrt(∑variance))
    resid = 0
    var = 0
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])] #see petrosian
        if (len(j) == 1 or len(j) == 2): continue
        j.append(i)
        L_ass = L[j]
        
        #determine rank
        L_rank = stat.rankdata(L_ass, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (len(L_ass))
        
        resid = resid + (rank - exp)
        var = var + (1/12.0 * (-1 + (len(L_ass)-1)**2))
        
        #troubleshooting
        if(i % 500 == 0): print i, resid, var
        
    t = resid / math.sqrt(var)
    return t

def tau_n(Z, L, Lmin): #petrosian said: normalize and see if there's a difference
    resid = 0
    var = 0
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] >= Lmin[i] and Z[m] <= Z[i])] #see petrosian
        rankIndex = np.where(np.array(j) == i)[0][0]
        L_ass = L[j]
        
        #determine rank
        L_rank = stat.rankdata(L_ass) #ranks of all luminosities
        rank = L_rank[rankIndex] #determine rank of data point i
        rank = rank / (1 + len(L_ass)) #normalize
        exp = 0.5
        
        resid = resid + (rank - exp)
        var = var + (1/12.0)
        
        #troubleshooting
        if(i % 500 == 0): print i, resid, var
        
    t = resid / math.sqrt(var)
    return t

def tau(Z, L, Lmin): #as defined in EP; tau = ∑(resid/var)
    t = 0
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] >= Lmin[i] and Z[m] <= Z[i])] #see petrosian
        rankIndex = np.where(np.array(j) == i)[0][0]
        L_ass = L[j]
        
        #determine rank
        L_rank = stat.rankdata(L_ass) #ranks of all luminosities
        rank = L_rank[rankIndex] #determine rank of data point i
        exp = 0.5 * (1 + len(L_ass))
        #print(len(L_ass))
        
        resid = (rank - exp)
        var = (1/12.0 * (len(L_ass)**2 - 1))
        
        if(len(L_ass) == 1):
            t = t
            #print 'smallest ass: ', i
        else:
            t = t + resid / math.sqrt(var)
        
        #troubleshooting
        if(i % 500 == 0): print i, resid, var
    
    t = t / math.sqrt(len(Z))
    return t

    
#tau(data_s[:,0], data_s[:,3], data_s[:,4])

# In[15]: Local Luminosity creation, given k:

#g(z) such that L' = L/g(z). Taken from Singal et al. (2013) and considers high-redshift objects.
def g(z, k):
    z_cr = 3.7
    return (1 + z)**k /(1 + ((1 + z) / z_cr)**k)

def localize(Z, L, Lmin, k):
    L_local = [L[p] / g(Z[p],k) for p in range(0, len(L))]
    Lmin_local = [Lmin[p] / g(Z[p],k) for p in range(0, len(L))]
    return np.array(L_local), np.array(Lmin_local)

#test with k = 3:
k = 3.5
L_local, Lmin_local = localize(data[:,0], data[:,3], data[:,4], k)
print 'L_local:', L_local, '\n'

#graph of L' vs. z
plt.figure(2)
plt.plot(data[:,0], np.log10(L_local),'.', markersize=1, label="my data", color='black')
plt.title("log(L') vs. z for k = 3")

Lmin_local_sorted = Lmin_local[np.argsort(data[:,0])] #used to graph the red line
plt.plot(zmin, np.log10(Lmin_local_sorted), markersize=3, label="my data", color='red')
axes = plt.gca()
axes.set_xlim([0,6])
axes.set_ylim([29,33])

#tau(data[:,0], L_local, Lmin_local)

# In[16]: Choose a sample of data: 
data_s = data[::10,:]

# In[17]: Iterate through k's to graph tau vs k:

def tauvsk(K, data):
    s = []
    for k in K:
        L_loc, Lmin_loc = localize(data[:,0], data[:,3], data[:,4], k)
        print('\nk = ' +  str(k))
        t = tau_(data[:,0], L_loc, Lmin_loc)
        print ('\ntau = ' + str(t))
        s.append(t)
    return np.array(s)

K = np.arange(3., 4., 0.1)
Tau = tauvsk(K, data_s)
# In[18]: Tau vs k plot:
plt.figure(3)
plt.title("tau vs k")
plt.plot(K, Tau)
axes = plt.gca()
axes.set_xlim([min(K), max(K)])
axes.set_ylim([-3,3])
plt.plot(np.arange(0, 10), np.zeros(10), color = 'black')
plt.plot(np.arange(0, 10), np.zeros(10) + 1, '--', color = 'black', linewidth=1)
plt.plot(np.arange(0, 10), np.zeros(10) - 1, '--', color = 'black', linewidth=1)
plt.show()