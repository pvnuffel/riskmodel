# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:50:43 2015

@author: pieter
"""

import numpy as np
import matplotlib.pylab as plt
import sys, glob
import scipy


#Nlist = scipy.array([1000,2000,3000,4000,5000,6000,7000,8000,9000,10000])
Nlist = scipy.array([1000,2000,4000,8000,16000,32000])
#np.savetxt('results/data/Jv_ms_e-6.out' , Jv_ms_list)
#np.savetxt('results/data/rho_ms_e-6.out' , rho_ms_list)
Jv_ms_eps6 = np.loadtxt('results/data/cluster/Jv_ms_eps-6.out' )
#rho_ms_e6 = np.loadtxt('results/data/rho_ms_e-6.out' )
Jv_ms_eps5 = np.loadtxt('results/data/cluster/Jv_ms_eps-5.out' )
#rho_ms_e6_cl = np.loadtxt('results/data/cluster/rho_ms_e-6.out' )
#Jv_ms_e5 = np.loadtxt('results/data/Jv_ms_e-5.out' )
#rho_ms_e5 = np.loadtxt('results/data/rho_ms.out' )

#plot_ms_rho = plt.plot(Nlist, rho_ms_list)
plot6 = plt.plot(Nlist, Jv_ms_eps6, label = '$\epsilon=10^{-6}$')
#plot6cl = plt.plot(Nlist, Jv_ms_e6_cl, label = '$\epsilon=10^{-6}$')
plot5 = plt.plot(Nlist, Jv_ms_eps5, label = '$\epsilon=10^{-5}$')

plt.xscale('log')
plt.yscale('log')

#Jv_norm = np.loadtxt('results/data/Jv.out')

#Jv_var = scipy.zeros(M)
#average= Jv_norm[-1]
#for m in range(1,M):
#    a= np.power((Jv_norm[m] - average),2) 
#    Jv_var[m]= np.sqrt(a)

#ax=plt.subplot(111)
#ax.set_xlim(1, 5000)
#ax.set_ylim(0, np.max(Jv_var))

#plot_var = plt.plot(Jv_var)
plt.xlabel('$N$', fontsize = 16)
plt.ylabel(r'$\frac{\mathbb{E}[ (\hat{\mathbf{\rho} } - \mathbf{\rho}_{FP} )^2 ] } {|| \mathbf{1}||^2}    $', fontsize = 28)
plt.ylabel(r'$\mathbb{E}[ (\hat{\mathbf{\rho} } - \mathbf{\rho}_{FP} )^2 ]   $', fontsize = 18)
plt.ylabel(r'$\mathbb{E}[ (\hat{\mathbf{J} } \cdot \mathbf{v}  - \mathbf{J_{FP}} \cdot \mathbf{v})^2 ]   $', fontsize = 18)
plt.legend([plot5, plot6], loc='best')
plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
#plt.savefig("results/plots/ms_Jveps-5eps-6.pdf")
plt.show()
