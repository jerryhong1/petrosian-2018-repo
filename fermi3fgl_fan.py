#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# blazars from Fan et al. (2016), obtained from the 3FGL telescope at Fermi.
# In[0]: read and parse data; 999 sources with known redshift
from __future__ import division
import csv
import numpy as np
import quasars as quas
from quasars import QuasarData, Band
import matplotlib.pyplot as plt

Z = []
Lr = []
Lo = []
Lx = []
Lg = []
file_begin = 67

with open('fan2016.tsv') as tsv:
    i = 0
    for line in csv.reader(tsv, delimiter="\t"):
        if (i >= file_begin):
            line = [s.replace(" ", "") for s in line]
            Z.append(line[1])
            Lr.append(line[3])
            Lo.append(line[4])
            Lx.append(line[5])
            Lg.append(line[6])
        i = i + 1

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
Lr = str2float(Lr)
Lo = str2float(Lo)
Lx = str2float(Lx)
Lg = str2float(Lg)
z_index = np.nonzero(Z)[0] #indeces where Z is known; will create a quasar object based on it

bands = ['r', 'o' , 'x', 'g'] #optical, radio, x-ray, gamma
nus= {'r': 1.4e9,
      'o': 4.68e14,
      'x': 2.416e17,
      'g': 2.416e23} # frequencies, all in Hz (note that optical is NOT 2500 A, but rather 6400
                     # uses r-band -- i could convert, if i wanted to...)

# In[3]: Plot
i = np.where((Z != 0.) & (Lo != 0.))[0]
plt.plot(Z[i], Lo[i], '.', markersize = 2)

# In[4]: K-correction for gamma: integrate over the gamma band and...    

# In[2]: Convert from nuLnu to Lnu
Lr = Lr / nus['r']
Lo = Lo / nus['o']
Lx = Lx / nus['x']
Lg = Lg / nus['g']

fmin = 3e-12 #for gamma, integrated over the whole band; see 

