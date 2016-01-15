# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""

import tables
import scipy
import sys
import numpy as np
from pylab import *
import os


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
    sigma =1
    D= 1./2
    Dt = 1.
    seed = 16
    # discretization parameters
   #h=2e-2 # kde mesh width 
    dx = 5e-2
#    dxlist= [1e-2,  5e-3, 1e-3, 5e-4] 
  #  dxlist= [1e-2, 1e-3]
    #dt = min(dxlist)**2/(2*D)
    dt=2e-3
    r= (2.*D*dt)/(dx)**2
    print'r = ' ,  r

   # r= (2.*D*dt)/(min(dxlist))**2
    if r>1: 
        print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
        sys.exit("Ending Script. Make sure that r<1") 
    xL = -1.7
    xR = 1.7
    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha])
    param=pde.FokkerPlanck.getDefaultParameters()
    param['Dt'] = Dt
    param['dt']=dt         
  #  norm_Jv_fp= scipy.zeros(len(dxlist))
 #   norm_rho_fp= scipy.zeros(len(dxlist))
   # rho_Dt_fp[i]= np.ndarray(shape=) #[scipy.zeros((len
    

    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho
 
    fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
    rho_Dt=fp_pde.u_Dt

#    v=scipy.zeros_like(grid)
# #   v[200]=-1
# #   v[300]=1
#    for j in range(len(grid)): 
#        v[j]= np.sin(j*2*np.pi/len(grid))

    #  plt.show()    
    # CREATING LINEAR SOLVER
    gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
    gmres_param['tol']=1e-7
    gmres_param['print']='long'
    gmres_param['builtin']=True
    linsolv = GMRES.GMRESLinearSolver(gmres_param)
    linsolv2 = GMRES.GMRESLinearSolver(gmres_param)

    # CREATING NEWTON SOLVER
    newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
    newt_param['rel_tol']=1e-7
    newt_param['abs_tol']=1e-7
    newt_param['print']='short'
    newt_param['max_iter']=3
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
    nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)

    # CREATING POINT
    psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
    psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolver()

    # POINT PARAMETERS
    point_param = Point.Point.getDefaultParameters()
    point_param['artificial']=[2]
    point_param['artificial_condition']=[probability.condition]
    
 #   pprev = p
    points = []
    p = Point.Point(fp_pde,nsolv,None, point_param)
    print "Solve for fixpoint"
    p.correct()
    print "Found! Add to list"
    points.append(p)
    
#    sigma =2
#    lambd = scipy.array([a,sigma,alpha,beta,zeta])
#    fp_pde2 = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
#    p2 = Point.Point(fp_pde2,nsolv2,None,point_param)
#    p2.correct()
#    points.append(p2)            
#      
#    continuer_param= Continuer.Continuer.getDefaultParameters()
#    branch = Continuer.Continuer(points, continuer_param)
#    branch.bcontinue(4)
#           
#
#    
    
        #rho_Dt_fp = fp_pde.u_Dt
        #plt.plot( grid, rho_Dt_fp)
     #   plt.plot( rho_Dt_fp)  
     #   norm_rho_fp[i] = norm(rho_Dt_fp)/sqrt(len(rho_Dt_fp))
       # norm_Jv_fp[i] = norm(Jv_pde)/sqrt(len(Jv_pde))
        