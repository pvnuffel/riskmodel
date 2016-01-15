# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 17:58:57 2015

@author: pieter
"""

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
from matplotlib.path import Path #to draw trriangles indicating order of convergence
import matplotlib.patches as patches #to draw trriangles indicating order of convergence

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


if __name__=="__main__":

    Mlist= [ 100, 200, 300, 400, 500]
    
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000]) #,128000])
  #  Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000, 128000, 256000])
    Nlist=scipy.array([1000])
    Nlist_inv = 1.0/Nlist
    
        #SDE       

    E_rho = scipy.zeros(((len(Mlist), 34)))    # [m x N_i x dx] = [5 x 8 x 3400]
    for M_factor in range(0,len(Mlist)):
        m= (M_factor+1)*100
        E_rho[M_factor] = np.loadtxt('test-E_rho_m%d.out' %m)
#    sq_E_rho = scipy.zeros((len(Mlist), len(Nlist))) 
#    for M_factor in range(0,5):
#        m= (M_factor+1)*20
#        sq_E_rho[M_factor] = np.loadtxt('data/8-12-sq_E_rho_m%d.out' %m)

#    sde_Jv_sq =scipy.zeros((len(Mlist), len(Nlist)))
#    for M_factor in range(0,len(Mlist)):
#        m= (M_factor+1)*20
#        sde_Jv_sq[M_factor] = np.loadtxt('data/fin-Jv_sq_m%d.out' %m)
    E_Jv = scipy.zeros(((len(Mlist), 34)))    # [m x N_i x dx] = [5 x 8 x 3400]
    for M_factor in range(0,len(Mlist)):
        m= (M_factor+1)*100
        E_Jv[M_factor] = np.loadtxt('test-E_Jv_m%d.out' %m)     
#    sq_E_Jv = scipy.zeros((len(Mlist), len(Nlist)))  
#    for M_factor in range(0,len(Mlist)):
#        m= (M_factor+1)*20
#        sq_E_Jv[M_factor] = np.loadtxt('data/fin-sq_E_Jv_m%d.out' %m) 
        
#   rho_Dt_sde = E_rho[-1][-1]
    Jv_sde = E_Jv[-1]  #dim n_
    bins = len(Jv_sde)
    
    
    def resize( original_vector, new_size):
        resize_factor = int (len(original_vector)/new_size) 
        print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"  
        new_vector = scipy.zeros(new_size)
        for i in range (0,new_size):
            bin_av = 0
            for j in range (0,  resize_factor ):
                bin_av = bin_av+ original_vector[i*resize_factor+j]
                new_vector[i] = bin_av/resize_factor
        return new_vector
  

    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0+0.08, box.y0, box.width*0.95, box.height]) 
    #r = scipy.zeros(10)  
    rlist= np.array([0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1])
   # rlist= [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.1]
    bias_sq_rlist = scipy.zeros(len(rlist))    
    for r_i in range(0, len(rlist)):
      #  bias_sq_rlist[r_inv-1] =  norm(resize( np.loadtxt('Jv_dx=e-2_r_inv=%d.out' %r_inv)  , bins) - Jv_sde )**2  Jv_dx=e-2_r=%d.out
        bias_sq_rlist[r_i] =  norm(resize( np.loadtxt('Jv_dx=e-2_r=%.3f.out' %rlist[r_i])  , bins) - Jv_sde )**2
    rlist =rlist/10000
      #  r[r_inv-1]= 1.0/r_inv
    plt.plot( rlist, bias_sq_rlist /bins, label ='$\Delta x = 10^{-2} $', linewidth=2  ) 
    plt.xlabel('$\Delta t$', fontsize = 14)
    plt.ylabel(r'$\left(\hat{\mathbf{Bias}}(\mathbf{\hat{Jv}} , \mathbf{Jv}_{FP} \right))^2$'  , fontsize = 15)
   # plt.xscale(log_flag)
   # plt.yscale(log_flag)
   # plt.ticklabel_format(style='sci',  scilimits=(0,10))
  #  plt.ticklabel_format(style='sci') 
    plt.ticklabel_format(useOffset=False, axis='y') 
    plt.savefig('plots/ftcs_checkup_linear_in_delta_t.pdf')
    plt.show()

    
 
    

