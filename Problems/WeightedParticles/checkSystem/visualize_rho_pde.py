# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""

import tables
import scipy
import sys
import numpy as np


from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../..")

import pde

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability
import Continuer.Continuer as Continuer 




if __name__=="__main__":
   # D = 1./2.
    sigma =1.
    D = sigma**2
    D=1.
    Dt= 6e-1
    seed = 16
    # discretization parameters
   #h=2e-2 # kde mesh width 
    dx = 1e-2
#    dxlist= [1e-2,  5e-3, 1e-3, 5e-4] 
  #  dxlist= [1e-2, 1e-3]
   # dt = min(dxlist)**2/(2*D)
    dt = 1e-5
    r= (2.*D*dt)/(dx)**2
    #r= (2.*D*dt)/(min(dxlist))**2
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
  #  norm_Jv_fp= scipy.zeros(len(dxlist))
  #  norm_rho_fp= scipy.zeros(len(dxlist))
   # rho_Dt_fp[i]= np.ndarray(shape=) #[scipy.zeros((len
    

    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho
    fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)

    v=scipy.zeros_like(grid)
 #   v[200]=-1
 #   v[300]=1
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))
        
    rho = fp_pde.u_Dt

   # Jv_pde = fp_pde.applyJacobian(v)    

#    for t in range(1,int(Dtot/Dt)): 
#        print t
#        plt.plot(grid, fp_pde.u_Dt)
#        plt.savefig('movieplots/plot_rho_t%d.pdf' %t)
    #  plt.show()    

        #rho_Dt_fp = fp_pde.u_Dt
        #plt.plot( grid, rho_Dt_fp)
     #   plt.plot( rho_Dt_fp)  
     #   norm_rho_fp[i] = norm(rho_Dt_fp)/sqrt(len(rho_Dt_fp))
       # norm_Jv_fp[i] = norm(Jv_pde)/sqrt(len(Jv_pde))
        