"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
#import tables
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
    dx = 5e-4
   # dt=1e-6
    dt = dx**2/(2*D)
    print Dt/dt, " timesteps"
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
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho

    v=scipy.zeros_like(grid)
    v[200]=-1
    v[300]=1
    
    
    #PDE
    t0 = time.time()
    print 'Start solving pde'
  #  fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
 #   Jv_pde = fp_pde.applyJacobian(v)    
#    rho_Dt_fp = fp_pde.u_Dt
 #   np.savetxt('results/data/rho_fp_dx_5e-5.out', rho_Dt_fp)
 #   np.savetxt('results/data/rho_Jv_dx_5e-5.out', Jv_pde)
            
    #SDE
    t1 = time.time()
    print "Simulation time for solving pde: " , t1-t0, " seconds"
    
    print 'Start solving sde'
    dx_fp= dx
    dx = 1e-2
    factor = int (dx/dx_fp)
    print "Discretisation for solving sde is ", factor, " times coarser than the discretisation for solving the pde"
    dt = 1e-4 # ti
    N=1000
    Nlarge= N
    Nsmall =N
        
    grid = scipy.arange(xL+dx/2.,xR,dx)
    
    if(len(grid)>N):
        print "The number of particles is smaller than the number of bins! Please increase N or increase dx"
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
  #  print "rho: ", rho     
    
    
#    v=scipy.zeros_like(grid)
#    v[200]=-1
#    v[300]=1
#    

    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=0
    lifting =  inv_transform.Sampler(sampler_param)
    
    h=2e-2 # kde mesh width     
    M=10  #number of monte carlo steps                
    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=h
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
    param['discr_factor']= factor
    param['Nlarge']=N
    param['Nsmall']=N    
    param['Dt'] = Dt 
    param['dt']=dt
   # fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param=param)           
    
   # rho_tot = scipy.zeros(len(rho)*factor)
    rho_var = scipy.zeros(len(rho))
    Jv_sde = scipy.zeros(len(v))   
    Jv_norm = scipy.zeros((M+1))
    rho_norm = scipy.zeros((M+1))
    rho_norm = scipy.zeros((M+1))
    Nlist = scipy.array([1000,2000,3000,4000,5000,6000,7000,8000,9000,10000])
    Nlist = scipy.array([800,1600,3200, 6400])
    nN= len(Nlist)
    rho_ms_list = scipy.zeros(nN)
    Jv_ms_list =  scipy.zeros(nN)
     

    for n in range(len(Nlist)):  
        N = Nlist[n]     
        param['Nlarge']=N
        param['Nsmall']=N    
        
        print 'run simulation with N = ', N, 'particles'
        rho_ms=0             
        Jv_ms=0    
        for m in range(1,M+1):    #M steps
            print "running with seed ", m 
            sampler_param = inv_transform.Sampler.getDefaultParameters()
            sampler_param['seed']=m
            lifting = inv_transform.Sampler(sampler_param)      #  dit duurt lang                   
            fp_sde=None
            fp_sde = particles.Particles(lifting,restriction,rho,grid,lambd, param=param)
            print "calculate .u_Dt"
            rho_Dt = fp_sde.u_Dt 
            rho_ms = rho_ms + (norm(rho_Dt))**2
        #    rho_tot = rho_tot +rho_fine
      #  rho_norm[m] = norm(rho_tot/m - rho_Dt_fp) #better to do this in post-processing-step?
        
      #  print rho_tot/m
            print "calculate JVs"
            #Jv_sde =fp_sde.applyJacobian(v) 
           # Jv_ms = Jv_ms + (norm( Jv_sde  -Jv_pde))**2
    
#        Jv_norm[m]= norm(Jv_sde/m - Jv_pde)
     #   Jv_sde_sum = Jv_sde_sum + Jv_sde
     #  diff =  Jv_sde - Jv_pde
     #   totaldiff = totaldiff + diff
#        averagediff = totaldiff/(m)
#        diff_norm[m] = norm(averagediff)
#        Jv_norm[m] = norm(Jv_sde_sum/(m))       
        rho_ms_list[n]  =rho_ms/(M)
       # Jv_ms_list[n] =Jv_ms/(M)
#Jv_ms = Jv_ms/M
#rho_av = rho_tot/M
#Jv_sde = Jv_sde/M

#np.savetxt('rho_ms_e-6.out' , rho_ms_list)
#np.savetxt('Jv_ms_e-6.out' , Jv_ms_list)

print "End of simulation"
now = time.time()
print "Simulation time for solving pde: " , t1-t0, " seconds"
print "Simulation time for solving sde: " , now-t1, " seconds"
print "Total simulation time " , now-t0, " seconds"