# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:02:34 2015

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

from scipy.linalg import norm
from scipy import sqrt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")

import pde
import histogram
import particles
import inv_transform_sampling as inv_transform

if __name__=="__main__":
    D = 1./2.
    Dt = 1e-2
    seed = 16
    # discretization parameters
   #h=2e-2 # kde mesh width 
    dx = 1e-5
   # dt=1e-6
    dt = dx**2/(2*D)
    print Dt/dt, " timesteps"
    if(Dt/dt <5):
        print "Choose more timesteps"
    r= (2.*D*dt)/(dx)**2
    if r>1: 
        print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
        sys.exit("Ending Script. Make sure that r<1") 
    xL = -1.7
    xR = 1.7
    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
    param=pde.FokkerPlanck.getDefaultParameters()
    param['Dt'] = Dt
    param['dt']=dt         

    grid = scipy.arange(xL+dx/2.,xR,dx)
    print len(grid), " discretisation steps"
    
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
 #   print "rho: ", rho

    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))

  #  v=scipy.zeros_like(grid)

#    v=scipy.zeros_like(grid)
 #   v[20]=1
 #   v[30]=-1
        
    #v=scipy.zeros_like(grid)
#    for j in range(2995,3005): 
#        v[j]= 1
#    print "perturbation vector: ",  v

   # v[300]=1
    
    #PDE
    t0 = time.time()
    print 'Start solving pde'
    
    
    
    eps_list_exponents=[4, 5,6,7]
    eps_list= [1e-4, 1e-5, 1e-6, 1e-7]
    
    
    
    for eps_i in range(len(eps_list)):
        
        param['eps']= eps_list[eps_i]
        fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
        Jv_pde = fp_pde.applyJacobian(v)    
        rho_Dt_fp = fp_pde.u_Dt
        np.savetxt('25-11-Rho_pde_dx_e-4_eps-%d.out' %eps_list_exponents[eps_i] , rho_Dt_fp)
        np.savetxt('25-11-Jv_pde_dx_e-4_eps-%d.out' %eps_list_exponents[eps_i], Jv_pde)
            
        
    t1 = time.time()
    print "Simulation time for solving pde: " , t1-t0, " seconds"