"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np
import matplotlib.pylab as plt
from numpy import *

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


def get_fixpoint( mu, sigma, alpha):
    xi=0
    C=  sigma**2/(2*alpha)
    if (3*C>1):
        print "Make sure that 3D/(2alpha) <1 to get non-zero equilibrium solutions"
    else:
        xi = sqrt(1-3*C)*(1+ 6*mu/(sigma**2)*C**2*(1-2*C)/(1-3*C)) 
    return xi


if __name__=="__main__":
    
    sigma=1.
    mu = 1
    alpha= 5.

    xi_abs= get_fixpoint( mu, sigma, alpha)
    print 'Analytical solution = ' , xi_abs
       
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000])
       
       
    xL = -1.7
    xR = 1.7
    dx = 1e-2
    grid = scipy.arange(xL+dx/2.,xR,dx)
    dt = 0.01
    
    Dt = 10000
    timesteps = scipy.arange(dt,Dt+dt,dt)
    xi1 = full(Dt/dt, xi_abs)
    xi2 = full( Dt/dt, -xi_abs)
    dt = 1e-2 # ti
    
    alphas = scipy.arange(0, 6, 0.11)
    fixpoints = scipy.zeros_like(alphas)
    for i in range (len(fixpoints)):
        fixpoints[i] = get_fixpoint(mu, sigma, alphas[i])
    
    plot_bif = plt.plot(alphas, fixpoints)
    plot_bif = plt.plot(alphas, -fixpoints)
    plt.ylabel(r'$x*$', fontsize = 18)
    plt.xlabel(r'$\alpha$', fontsize = 16)
  #  plt.savefig('plots/bif_alpha.pdf')
    plt.show()
    
    
    rho_ss_alpha1= np.loadtxt('pfix_alpha1.txt')
    xmean_alpha1 = np.dot(grid,rho_ss_alpha1)*dx
    rho_ss_alpha5 =np.loadtxt('pfix_alpha5.txt')
    xmean_alpha5 = np.dot(grid,rho_ss_alpha5)*dx
    plt.axvline(xmean_alpha1, linewidth=2)
    plot_anal =plt.axvline( get_fixpoint( mu, sigma, 1), color='red',  linestyle='--', linewidth=2, label =  r'${\bar{x}}_{analytical}$')

    plot_stable = plt.plot(grid, rho_ss_alpha1, label = r'${\rho^*}_{stable}$')
    plt.axvline(xmean_alpha5, color='green', linewidth=2)
    plt.axvline( -get_fixpoint( mu, sigma, 5), color='red',  linestyle='--', linewidth=2)
    plt.xlabel('$x$', fontsize = 16)
    plt.ylabel(r'$\rho*$', fontsize = 18)
    delta_x_alpha1 = np.abs(xmean_alpha1- get_fixpoint( mu, sigma, 1))
    delta_x_alpha5  =  np.abs(xmean_alpha5 + xi_abs)
    plt.title( r'$\Delta \bar{x} (\alpha=5) =%.3f}$    ' %delta_x_alpha5 + r'$\Delta \bar{x} (\alpha=1) =%.3f}$' %delta_x_alpha1 , fontsize=11)
    plot_meta= plt.plot(grid, rho_ss_alpha5, label = r'${\rho^*}_{metastable}$')
    plt.legend([plot_anal, plot_meta, plot_stable], loc='best')

   # plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
    plt.legend(loc=0, prop={'size':11}) #fontsize = 'x-small')
    plt.gca().set_ylim([0, 1.5])   
 #   plt.savefig('plots/steady_states.pdf')
    plt.show()    
    
    
    print 'Analytical solution: systemic risk = ', xi_abs
    print 'Newton-Krylov solution: systemic risk = ', xmean_alpha5
  
    
   # plot_sysrisk =  plt.plot( timesteps, xmean  , label = '$\$')
    plot_sysrisk = plt.plot( timesteps, xi1, linestyle='--', linewidth = 2 , color = 'red', label = '$\$')
    plot_sysrisk = plt.plot( timesteps, xi2, linestyle='--', linewidth = 2,  color = 'red' , label = '$\$')
    #plot_rho =  plt.plot(rho_Dt_sde, label = '$\$')
    plt.ylabel(r' $\bar{x}$', fontsize = 18)
    plt.xlabel('$t$', fontsize = 16)
    plt.title(r'$N=100$, $\mu=%.1f$, $\sigma=%d$, $\alpha=%d$' %(mu, sigma, alpha) )
 #  plt.savefig('plots/xmean_alpha%d.pdf'%alpha, bbox_inches='tight')
    plt.show()
       
