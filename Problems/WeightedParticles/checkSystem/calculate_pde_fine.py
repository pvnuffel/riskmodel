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

import matplotlib.pylab as plt

if __name__=="__main__":
    D = 1./2.
    Dt = 1e-2
    seed = 16
    # discretization parameters
   #h=2e-2 # kde mesh width 
    dx = 1e-2
    dt=1e-6
#    dt = dx**2/(2*D)*r
    print Dt/dt, " timesteps"
#    r= (2.*D*dt)/(dx)**2
#    print 'r= ', r
#    if r>1: 
#         print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
#         sys.exit("Ending Script. Make sure that r<1") 
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
 
    #PDE
    t0 = time.time()
    print 'Start solving pde'
    
    
    #rlist= [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]
    dxlist= [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
    for dx_i in range (0,len(dxlist)):
        dx =dxlist[dx_i]
        print 'dx= ', dx
        r= (2.*D*dt)/(dx)**2
        print 'r= ', r
        if r>1: 
            print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
            #sys.exit("Ending Script. Make sure that r<1")       
            
        grid = scipy.arange(xL+dx/2.,xR,dx)
        print len(grid), " discretisation steps"
        rho = scipy.ones_like(grid)
        rho = rho/(sum(rho)*dx)    
        
        v=scipy.zeros_like(grid)
        for j in range(len(grid)): 
            v[j]= np.sin(j*2*np.pi/len(grid))
        plt.plot(v)
    
        fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
        Jv_pde = fp_pde.applyJacobian(v)    
        rho_Dt_fp = fp_pde.u_Dt
        np.savetxt('rho_dx=%.3f.out' %dx, rho_Dt_fp)
        np.savetxt('Jv_dx=%.3f.out'%dx, Jv_pde)
            
    t1 = time.time()
    print "Simulation time for solving pde: " , t1-t0, " seconds"