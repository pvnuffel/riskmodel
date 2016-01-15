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

import pde
import histogram
import particles
import inv_transform_sampling as inv_transform

if __name__=="__main__":
    D = 1./2.
    Dt = 1e-2
    # discretization parameters
   #h=2e-2 # kde mesh width 
    xL = -1.7
    xR = 1.7
    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
  
    
    #PDE
#    print 'Start solving pde'
#    fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
#    Jv_pde = fp_pde.applyJacobian(v)    
#    rho_Dt_fp = fp_pde.u_Dt
 #   Jv_pde= np.loadtxt('Jv_fp_dx_1e-4.out')
#    rho_Dt_fp= np.loadtxt('rho_fp_dx_1e-4.out') 
 #   np.savetxt('rho_fp_dx_5e-5.out', rho_Dt_fp)
 #   np.savetxt('rho_Jv_dx_5e-5.out', Jv_pde)
            
    #SDE
    t1 = time.time()
    
    print 'Start solving sde'
    dx = 1e-1

    dt = 1e-4 # ti
    N=1000
    Nlarge= N
    Nsmall =N
    
    
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    
    

    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))

#    v=scipy.zeros_like(grid)
#    v[20]=1
#    v[30]=-1
     
    
  #  sampler_param = inv_transform.Sampler.getDefaultParameters()
  #  sampler_param['seed']=0
  #  lifting =  inv_transform.Sampler(sampler_param)
    
    h=2e-2 # kde mesh width     
              
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

   # Nlist = scipy.array([1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000])
    Nlist = scipy.array([1000]) #, 800, 1600, 3200])
    
    param['eps']= 1e-5 
    M=100 #number of monte carlo steps  
    
    nN= len(Nlist)

    rho_sq_sum = scipy.zeros(nN)
    Jv_sq_sum = scipy.zeros(nN)
    E_rho_sum= scipy.zeros((nN,len(rho)))
    E_Jv_sum= scipy.zeros((nN,len(rho)))
    
    #Jv_sq_sum =  np.loadtxt('test-Jv_sq_m100.out')*100    #use realizations from previous simulations
  #  E_Jv_sum= np.loadtxt('fin-E_Jv_m100.out')*100
        
    
    for M_factor in range(0,5):     #iterate 5 times but remember MC-data
        
        for m in range(1,M+1):    #M steps
            print "running with seed ", M_factor*M  +m
            sampler_param = inv_transform.Sampler.getDefaultParameters()
            sampler_param['seed']=  M_factor*M  +m
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
                #rho_fine = fp_sde.make_rho_fine(rho_Dt) 
                rho_sq_sum[n] = rho_sq_sum[n] + (norm(rho_Dt))**2
                E_rho_sum[n] = E_rho_sum[n] + rho_Dt           #matrix
        
                print "calculate JVs"
                Jv =fp_sde.applyJacobian(v)                                                 
                Jv_sq_sum[n] = Jv_sq_sum[n] + (norm(Jv))**2
                E_Jv_sum[n] = E_Jv_sum[n] + Jv 
                
                x_prev_sim = fp_sde.x_Dt
              #  print 'len x-vector= ' , len( fp_sde.x_Dt) 
                w_prev_sim = fp_sde.w_prev
                #print w_prev_sim
                                 
                                 
        rho_sq = scipy.zeros(nN)
        Jv_sq = scipy.zeros(nN)
        sq_E_rho=scipy.zeros(nN)
        sq_E_Jv=scipy.zeros(nN)
        E_rho= scipy.zeros((nN,len(rho)))
        E_Jv= scipy.zeros((nN,len(rho)))        
        Mtot= M*(M_factor+1)
              
        for n in range(len(Nlist)): 
            rho_sq[n]=  rho_sq_sum[n]/Mtot
            E_rho[n] = E_rho_sum[n]/Mtot
            sq_E_rho[n]= (norm(E_rho[n]))**2
            
            Jv_sq[n]=  Jv_sq_sum[n]/Mtot
            E_Jv[n] = E_Jv_sum[n]/Mtot
            sq_E_Jv[n]= (norm(E_Jv[n]))**2
            
       # write_to_path = '/home/pieter/fokkerplanck/Problems/WeightedParticles/checkSystem/data/newfile_m%d.out' %Mtot
        #file1 = open('data/newfile_m100.out', "w")
        np.savetxt('test-E_Jv_m%d.out'%Mtot, E_Jv)
        np.savetxt('test-E_rho_m%d.out'%Mtot, E_rho)
       #np.savetxt('5-12-E_rho_m%d.out' %Mtot, E_rho )
       # np.savetxt('5-12-sq_E_rho_m%d.out' %Mtot,  sq_E_rho )        
       # np.savetxt(file1,  '42')
       # file1.close()
       #np.savetxt('test-E_Jv_m%d.out' %Mtot, E_Jv )
      #  np.savetxt('test-sq_E_Jv_m%d.out' %Mtot, sq_E_Jv)
    

    print "End of simulation"
    now = time.time()
  #  print "Simulation time for solving pde: " , t1-t0, " seconds"
    print "Simulation time for solving sde: " , now-t1, " seconds"
#    print "Total simulation time " , now-t0, " seconds"