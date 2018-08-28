#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from quasars import Band, QuasarData
import quasars as quas

#quasarPlotter

# makes everything but the title and saving it
sigma = [-1, 1]
def tauvsk(band):
    b = band
    k = b.k_g
    K = b.k_array
    Tau = b.tau_array
    
    plt.figure()
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
    
    
def rvsalpha(Alpha, R):
    #alpha vs r
    plt.figure()
    plt.plot(Alpha, R, color = 'black', linewidth = 1)
    #plt.plot(Alpha[::10], R[::10], marker = '.', markersize = 8)
    alpha = np.interp(0, -np.array(R), Alpha)
    plt.plot([alpha, alpha], [-4, 0], color = 'red', linewidth = 1)
    plt.ylabel(r"Correlation Coefficient")
    plt.xlabel(r"$\alpha$")
    plt.text(alpha + 0.02, 0.01, r"$\alpha = $ " + str(round(alpha, 2)), color = 'red', fontsize = 20)
    
    axes = plt.gca()
    plt.plot(np.arange(-1, 3), np.zeros(4), color = 'black', linewidth = 1)
    axes.set_xlim([min(Alpha), 1])
    axes.set_ylim([min(R), max(R)])
