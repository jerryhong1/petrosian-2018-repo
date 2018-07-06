#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 10:37:24 2018

@author: jerryhong
"""

#debugging tau

import numpy as np
import scipy as sp
from scipy.integrate import quad
import scipy.stats as stat
import math
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits

data_s = data[::2000, :]

#L vs z
zmin = data_s[:,0][np.argsort(data_s[:,0])]
Lmin_sorted = data_s[:,4][np.argsort(data_s[:,0])]
        
#plt.plot(np.arange(0.01,5,0.05), np.log10(Lmax),'-', markersize=3, label="max values", color='red')


def tau_(Z, L, Lmin): #as defined in singal, petrosian papers in 2010s. tau = (∑resid)/(sqrt(∑variance))
    resid = 0
    var = 0
    for i in range(0, len(Z)):
        #create associated sets
        j = [m for m in range(0, len(Z)) if (L[m] > Lmin[i] and Z[m] < Z[i])] #see petrosian

        j.append(i)
        L_ass = L[j]
        Z_ass = Z[j]
        
        #associated set (see above cell)
        plt.figure()
        plt.plot(data_s[:,0], np.log10(data_s[:,3]),'.', markersize=1, label="my data", color='black')
        #Lmin vs z
        #plt.plot(zmin, np.log10(Lmin_sorted),'-', markersize=3, label="min values", color='red')


        plt.plot(Z_ass, np.log10(L_ass),'.', markersize=1, label="associated set for source", color='red')
        plt.plot(Z_ass[-1], np.log10(L_ass[-1]), '.', markersize=15, label="source", color = 'green')
        plt.plot(Z[i], np.log10(L[i]), '.', markersize=9, label="source", color = 'blue')
        plt.show()
        print(len(j))
        if(raw_input('continue?') == 'n'): break
        
        #determine rank
        if (len(j) == 1 or len(j) == 2): continue
        L_rank = stat.rankdata(L_ass, 'max') #ranks of all luminosities
        rank = L_rank[-1] - 1 #determine rank of data point i
        exp = 0.5 * (len(L_ass))
        
        resid = resid + (rank - exp)
        var = var + (1/12.0 * (-1 + (len(L_ass)-1)**2))
        
        #troubleshooting
        print i, rank, resid, var, resid/math.sqrt(var)
        
    t = resid / math.sqrt(var)
    return t

tau_(data_s[:,0], data_s[:,3], data_s[:,4])