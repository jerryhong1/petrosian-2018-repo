# coding: utf-8

# compare the USNO and SDSS v-band data provided by 3LAC 
# and the SDSS i-band data (SDSS_missing) provided by NED with each other.
# In[0]: read and parse data; 999 sources with known redshift
from __future__ import division
import csv
import numpy as np
import quasars as quas
from quasars import QuasarData, Band
import matplotlib.pyplot as plt
plt.style.use('mystyle.mplstyle')

names = [] # 3FGL names
cnames = [] #companion name
Z = [] # redshifts
Fr = [] # radio fluxes (mJy)
Fx = [] # x-ray fluxes (erg/s/cm^2)
V_USNO = [] # optical V magnitudes
V_SDSS = [] # optical V magnitudes
I_missing = [] # optical I magnitudes to be filled in
G_missing = []
R_missing = []

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
                Fx.append(line[11])
                V_USNO.append(line[12])
                V_SDSS.append(line[13])
                I_missing.append(0.0)
                G_missing.append(0.0)
                R_missing.append(0.0)
                #find magnitudes on NED for objects without optical data provided
#                if(line[12].strip() and line[13].strip()): 
#                print line[2]
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
            R_missing[index] = float(linesplit[11])
            G_missing[index] = float(linesplit[9])
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
R_missing = np.array(R_missing)
G_missing = np.array(G_missing)

# In[2]: Optical data analysis.
f0_i = 3.631e-20  # from SDSS paper (2010)
f0_v = 3.8363e-20  # works for USNO?

def magtoflux(V, f0):
    F = []
    for v in V:
        if v == 0.0:
            F.append(0.0)
        else:
            F.append(quas.magtoflux(v, f0))
    return F
    
F_USNO = magtoflux(V_USNO, f0_v)
F_SDSS = magtoflux(V_SDSS, f0_v)
F_SDSS_i = magtoflux(I_missing, f0_i)

# assume V band is 550 nm: convert everything to optical
F_USNO = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_USNO]
F_SDSS = [quas.bandtoband(f, quas.lambda_v, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS]
F_SDSS_i = [quas.bandtoband(f, quas.lambda_i, quas.lambda_opt, quas.alpha_opt) for f in F_SDSS_i]

fmin_USNO = quas.bandtoband(quas.magtoflux(21, f0_v), quas.lambda_v, quas.lambda_opt, quas.alpha_opt)
fmin_SDSS = quas.bandtoband(0.08317e-26, quas.lambda_i, quas.lambda_opt, quas.alpha_opt)

# adopt Richards et al. (2006) k-correction for both
o_USNO = Band('o', fmin_USNO, F_USNO, quas.k_opt)
o_SDSS = Band('o', fmin_SDSS, F_SDSS, quas.k_opt)
o_SDSS_i = Band('o', fmin_SDSS, F_SDSS_i, quas.k_opt)

USNO = QuasarData(Z, [o_USNO])
SDSS = QuasarData(Z, [o_SDSS])
SDSS_i = QuasarData(Z, [o_SDSS_i])
USNO.sort()
SDSS.sort()
SDSS_i.sort()

# In[3]: L vs. z Plots
plt.figure(figsize=(8, 6))
plt.semilogy(USNO.Z, o_USNO.L, '.', markersize = 2)
plt.semilogy(USNO.Z, o_USNO.Lmin)
plt.title(r"3LAC USNO optical data + truncation ($V < 21)$")
plt.xlabel("z")
plt.ylabel("log(L)")

plt.figure(figsize=(8, 6))
plt.semilogy(SDSS.Z, o_SDSS.L, '.', markersize = 2)
plt.semilogy(SDSS.Z, o_SDSS.Lmin)
plt.title(r"3LAC SDSS optical data + truncation ($i < 19.1$)")
plt.xlabel("z")
plt.ylabel("log(L)")

plt.figure(figsize=(8, 6))
plt.semilogy(SDSS_i.Z, o_SDSS_i.L, '.', markersize = 2)
plt.semilogy(SDSS_i.Z, o_SDSS_i.Lmin)
plt.title(r"NED SDSS i-band optical data + truncation ($i < 19.1$)")
plt.xlabel("z")
plt.ylabel("log(L)")

# In[4]: Band to band comparison

# USNO vs. SDSS
both = np.where((V_USNO != 0) & (V_SDSS != 0))
plt.figure(figsize=(8, 6))
plt.plot(V_USNO[both], V_SDSS[both], '.')
plt.plot(range(12, 23), range(12, 23))
plt.title("V band observation, USNO vs. SDSS")
plt.xlabel("USNO")
plt.ylabel("SDSS")
axes = plt.gca()
axes.set_xlim([13,22])
axes.set_ylim([13,22])


both = np.where((o_USNO.F != 0) & (o_SDSS.F != 0))
plt.figure(figsize=(8, 6))
plt.plot(np.log10(o_USNO.F[both]), np.log10(o_SDSS.F[both]), '.')
plt.plot(range(-33, 3), range(-33, 3))
plt.title("Optical flux, USNO vs. SDSS provided in 3LAC")
plt.xlabel("USNO")
plt.ylabel("SDSS")
axes = plt.gca()
axes.set_xlim([-28, -24])
axes.set_ylim([-28, -24])

# SDSS 3LAC vs. SDSS NED
both = np.where((o_SDSS_i.F != 0) & (o_SDSS.F != 0))
plt.figure(figsize=(8, 6))
plt.plot(np.log10(o_SDSS_i.F[both]), np.log10(o_SDSS.F[both]), '.')
plt.plot(range(-33, 3), range(-33, 3))
plt.title("Optical flux, SDSS from NED (i band) vs. SDSS provided in 3LAC (v band)")
plt.xlabel("SDSS NED")
plt.ylabel("SDSS 3LAC")
axes = plt.gca()
axes.set_xlim([-28, -24])
axes.set_ylim([-28, -24])


both = np.where((V_SDSS != 0) & (G_missing != 0))
plt.figure(figsize=(8, 6))
plt.plot(V_SDSS[both], G_missing[both], '.')
plt.plot(range(12, 23), range(12, 23))
plt.title("V vs. g band observation, SDSS")
plt.xlabel("V")
plt.ylabel("g")
axes = plt.gca()
axes.set_xlim([12,22])
axes.set_ylim([12,22])

plt.figure(figsize=(8, 6))
plt.plot(V_SDSS[both], R_missing[both], '.')
plt.plot(range(12, 23), range(12, 23))
plt.title("V vs. r band observation, SDSS")
plt.xlabel("V")
plt.ylabel("r")
axes = plt.gca()
axes.set_xlim([12,22])
axes.set_ylim([12,22])

plt.figure()
plt.plot(R_missing, G_missing, '.')
plt.plot(range(12, 23), range(12, 23))
axes = plt.gca()
axes.set_xlim([12,22])
axes.set_ylim([12,22])