import sys

sys.path.append('../')
sys.path.append('../../')
sys.path.append('../../../')
sys.path.append('../system/')
sys.path.append('../lifting/')
sys.path.append('../restriction/')
sys.path.append('../modules/')

import scipy
import scipy.linalg
import numpy as np

import time 
import tables

import inv_transform
#import kde
import histogram
import particles
import pde
import precond

# Fokker Planck for SDE
a = 1  #advection coefficient
D =  1./2. #difussion coefficient
lambd = scipy.array([a,D])
Dt = 0.01 
N = 10000
eps = 1e-5

# discretization parameters
#h=2e-2 # kde mesh width 
dx = 0.1
dt = 1e-4

xL = -1.7
xR = 1.7
grid = scipy.arange(xL+dx/2.,xR,dx)
rho = scipy.ones_like(grid)
rho = rho/(sum(rho)*dx)
#print ("rho: ", rho)


v=scipy.zeros_like(grid)
for j in range(len(grid)): 
    v[j]= np.sin(j*2*np.pi/len(grid))

t0 = time.time()
print 'Start solving pde'

## Fokker Planck for SDE
#sampler_param = inv_transform.Sampler.getDefaultParameters()
#sampler_param['seed']=0
#lifting = inv_transform.Sampler(sampler_param)

restriction = histogram.Histogram()

param=particles.Particles.getDefaultParameters()
param['Nlarge']=N
param['Nsmall']=N    
param['Dt'] = Dt 
param['dt']=dt

param['eps']= eps
M=500 #number of monte carlo steps  

E_rho =scipy.zeros_like(rho)
E_Jv =scipy.zeros_like(rho)

for m in range(1,M+1): 
    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=m
    lifting = inv_transform.Sampler(sampler_param)   
    fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param)      
    # testing the run
    print (" run  ----")
    rho_Dt = fp_sde.integrate(rho,lambd)
    E_rho = E_rho + rho_Dt
    print (" perturbed run ---")
   # rho_Dt_pert = fp_sde.integrate(rho+eps*v,lambd)
    print ("Jacobian --- ")
    Jv_sde = fp_sde.applyJacobian(v)
    E_Jv = E_Jv + Jv_sde 

    # # testing the sampler algorithm :
    # 
    # rho_Dt = scipy.zeros_like(grid)
    # rho_Dt[5:-5] = 1.
    # rho_Dt = rho_Dt/sum(rho_Dt)*(grid[1]-grid[0])
    # xx = sampler.sample_density(rho_Dt,grid,N)
    # rho_Dt_resampled = fp_sde.restrict_histogram(xx)
    

np.savetxt('E_rho_M100_nw.out' , E_rho/M )
np.savetxt('E_Jv_M100_nw.out' , E_Jv/M )

print "End of simulation"
now = time.time()
  #  print "Simulation time for solving pde: " , t1-t0, " seconds"
print "Simulation time for solving sde: " , now-t0, " seconds"

