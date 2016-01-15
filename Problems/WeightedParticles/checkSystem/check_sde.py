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

import matplotlib.pylab as plt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")

import histogram
import particles
import inv_transform_sampling as inv_transform

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability
import Continuer.Continuer as Continuer 


if __name__=="__main__":
    sigma = 1.
    Dt = 100
    seed = 42
    # discretization parameters
   #h=2e-2 # kde mesh width 
   # dt=1e-6
    xL = -1.7
    xR = 1.7
    
    a = 0.1
    zeta = 0.
    alpha= 5
    beta = -1
    lambd = scipy.array([a,sigma,alpha,beta,zeta])
            
    #SDE
    t1 = time.time()
    
    print 'Start solving sde'
    dx = 1e-2
#    factor = int (dx/dx_fp)
#    print "Discretisation for solving sde is ", factor, " times coarser than the discretisation for solving the pde"
    dt = 1e-2 # ti
    N=10000
    Nlarge= N
    Nsmall =N
    
    
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho_left = scipy.ones(len(grid)/2)
    rho_right= scipy.ones(len(grid)/2)/5
#    rho_right= scipy.zeros(len(grid)/2)
    rho =scipy.r_[rho_left, rho_right] 
    rho = rho/(sum(rho)*dx)
    print rho

    
    

    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))

    
  #  sampler_param = inv_transform.Sampler.getDefaultParameters()
  #  sampler_param['seed']=0
  #  lifting =  inv_transform.Sampler(sampler_param)
    
    h=2e-2 # kde mesh width     
    M=1 #number of monte carlo steps                
    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=h
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
    param['Nlarge']=N
    param['Nsmall']=N    
    param['Dt'] = Dt 
    param['dt']=dt
#    param['eps']=1e-5
   # fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param=param)       

  #  Nlist = scipy.array([1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000])
    Nlist = scipy.array([N]) #, 800, 1600, 3200])
    eps_list_exponents=[5]
    eps_list= [1e-5]
    param['eps']=eps_list[-1]
   # Nlist = scipy.array([5,10, 200]) #,400])
  #  Nlist = scipy.array([4000])   
    
    nN= len(Nlist)

    rho_sq = scipy.zeros(nN)
    Jv_sq = scipy.zeros(nN)
    sq_E_rho=scipy.zeros(nN)
    sq_E_Jv=scipy.zeros(nN)
    E_rho= scipy.zeros((nN,len(rho)))
    E_Jv= scipy.zeros((nN,len(rho)))
    points = []

    for m in range(1,M+1):    #M steps
        print "running with seed ", m 
        sampler_param = inv_transform.Sampler.getDefaultParameters()
        sampler_param['seed']=m
        lifting = inv_transform.Sampler(sampler_param)         
        #fp_sde=None
        x_prev_sim= None
        w_prev_sim = None
        rg_state = None
        
        for n in range(nN):  
            N = Nlist[n] - np.sign(n)*Nlist[n-1]   #reduce the number of particle simulations by using the data from previous simulations
            if(len(grid)>N):
                print "The number of particles is smaller than the number of bins! Please increase N or increase dx"

            param['Nlarge']=Nlist[n]
            param['Nsmall']=N    
           # param['Nsmall']=Nlist[n]    
        
            print 'run simulation with N = ', N, ' = ', Nlist[n], ' - ' , np.sign(n)*Nlist[n-1]  , 'particles'
            if(rg_state != None): 
                lifting.set_rg_state(rg_state)
            fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, x_prev_sim, w_prev_sim, param=param )   
            rg_state = lifting.get_rg_state()  #We remember the internal state of the random generator to get new random number for sampling the new particles
            
            print "calculate .u_Dt"
            rho_Dt = fp_sde.u_Dt 
            xmean = fp_sde.x_mean
            #rho_fine = fp_sde.make_rho_fine(rho_Dt) 
         #   rho_sq[n] = rho_sq[n] + (norm(rho_Dt))**2
            E_rho[n] = E_rho[n] + rho_Dt           #matrix
            
            t0=time.time()
            print "Simulation time for solving sde: " , t0-t1, " seconds"
    
            print "calculate JVs"
         #  Jv =fp_sde.applyJacobian(v)                                                 
           # Jv_sq[n] = Jv_sq[n] + (norm(Jv))**2
            #E_Jv[n] = E_Jv[n] + Jv 
            

            #print w_prev_sim
            
                
                        # CREATING rho = scipy.ones_like(gridLINEAR SOLVER
            gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
            gmres_param['tol']=1e-7
            gmres_param['print']='short'
            gmres_param['builtin']=True
            linsolv = GMRES.GMRESLinearSolver(gmres_param)
            linsolv2 = GMRES.GMRESLinearSolver(gmres_param)
    
            # CREATING NEWTON SOLVER
            newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
            newt_param['rel_tol']=1e-7
            newt_param['abs_tol']=1e-7
            newt_param['print']='short'
            newt_param['max_iter']=4
            nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
            nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)
    
            # CREATING POINT
            psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
            psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolver()
    
            # POINT PARAMETERS
            point_param = Point.Point.getDefaultParameters()
            point_param['artificial']=[4]
            point_param['artificial_condition']=[probability.condition]
            
         #   pprev = p
            p = Point.Point(fp_sde,nsolv,psolver_im, point_param)
            p.correct()
            points.append(p)

            
            alpha=4
            lambd = scipy.array([a,sigma,alpha,beta,zeta])
            fp_sde2 = particles.Particles(lifting,restriction,rho,grid, lambd, x_prev_sim, w_prev_sim, param=param )   
            p2 = Point.Point(fp_sde2,nsolv2,psolver_im2,point_param)
            p2.correct()
            points.append(p2)            
            
            x_prev_sim = fp_sde.x_Dt
            w_prev_sim = fp_sde.w_prev
   
    for n in range(len(Nlist)): 
       # rho_sq[n]=  rho_sq[n]/M
        E_rho[n] = E_rho[n]/M
      #  sq_E_rho[n]= (norm(E_rho[n]))**2
        
      #  Jv_sq[n]=  Jv_sq[n]/M
      #  E_Jv[n] = E_Jv[n]/M
      #  sq_E_Jv[n]= (norm(E_Jv[n]))**2
        
        
    continuer_param= Continuer.Continuer.getDefaultParameters()
    branch = Continuer.Continuer(points, continuer_param)
    #branch.bcontinue(4)
       
        
    

    print "End of simulation"
    now = time.time()
    rho = rho/(sum(rho)*dx)
  #  print "Simulation time for solving pde: " , t1-t0, " seconds"

    print "Total simulation time " , now-t0, " seconds"