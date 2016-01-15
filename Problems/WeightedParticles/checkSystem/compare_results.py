"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np

from scipy.linalg import norm
from scipy import sqrt

#import matplotlib.pylab as plt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")


if __name__=="__main__":
    D = 1./2.
    Dt = 1e-2    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
    param['eps'] = 1e-5
    param['Dt'] = Dt
    param['dt']=dt         


    Jv_pde= np.loadtxt('Jv_fp_dx_1e-4.out')
    rho_Dt_fp= np.loadtxt('rho_fp_dx_1e-4.out') 


    
   # rho_sde = np.loadtxt('Jv_fp_dx_1e-4.out')       
    #Jv_sde = np.loadtxt('Jv_fp_dx_1e-4.out')
            
    #SDE
    t1 = time.time()
    print "Simulation time for solving pde: " , t1-t0, " seconds"
    
    print 'Start solving sde'
    dx_fp= dx
    dx = 1e-2
    rho_Dt_sde = rho_Dt
   # factor = int (dx/dx_fp)
    
    resize_factor = int (len(rho_Dt_fp)/len(rho_Dt_sde))
    print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"
    
    rho_coarse = scipy.zeros(len(rho_Dt_sde))
    for i in range (0,len(rho_coarse)):
        bin_av = 0
        for j in range (0,  resize_factor ):
            bin_av = bin_av+ rho_Dt_fp[i*resize_factor+j]
        rho_coarse[i] = bin_av/resize_factor
        
    pde_norm = norm(rho_coarse)**2
    print ( rho_tot  - pde_norm)
    
               
                        
        
  
    #dt = 1e-4 # ti
 
    
#    v=scipy.zeros_like(grid)
#    v[200]=-1
#    v[300]=1
#    
#
#    sampler_param = inv_transform.Sampler.getDefaultParameters()
#    sampler_param['seed']=0
#    lifting =  inv_transform.Sampler(sampler_param)
#    
#    h=2e-2 # kde mesh width     
#    M=10  #number of monte carlo steps                
#    param_histogram = histogram.Histogram.getDefaultParameters()
#    param_histogram['h']=h
#    restriction = histogram.Histogram(param_histogram)
#
#    param=particles.Particles.getDefaultParameters()
#    param['discr_factor']= factor
#    param['Nlarge']=N
#    param['Nsmall']=N    
#    param['Dt'] = Dt 

    
