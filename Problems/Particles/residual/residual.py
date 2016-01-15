import sys

sys.path.append('../')
sys.path.append('../../')
sys.path.append('../../../')
sys.path.append('../system/')
sys.path.append('../lifting/')
sys.path.append('../restriction/')

import scipy
import scipy.linalg

import inv_transform
#import kde 
import histogram 

import particles
import pde
import precond

# system parameters
a = 1
D = 0.1
lambd = scipy.array([a,D])
xL = -2.
xR = 2.

# method parameters
Dt = 0.05
eps_jac = 1e-3
dx = 0.05
grid = scipy.arange(xL+dx/2.,xR,dx)

# initial condition
rho0 = scipy.ones_like(grid)/(xR-xL)
rho = scipy.array(rho0)


# PDE-based system
param=pde.FokkerPlanck.getDefaultParameters()
print "param", param
param['D']=D
param['a']=a
param['Dt'] = Dt
param['eps_jac']=eps_jac
fp_pde = pde.FokkerPlanck(rho0,grid,pde.doublewell,lambd,param)
rho_Dt = fp_pde.u_Dt

# particle-based system

#param_kde = kde.KDE.getDefaultParameters()
#param_kde['h']=2e-2
#restriction = kde.KDE(param_kde)
restriction = histogram.Histogram()

param=particles.Particles.getDefaultParameters()
param['Dt'] = Dt
precond_param=precond.Precond.getDefaultParameters()
precond_param['Dstar']=0.
precond_param['sigma']=1.
param['precond']=precond.Precond(precond_param)

Nlist = scipy.array([10000,20000,50000,100000,200000,500000,700000])
nN = len(Nlist)
tolerance0 = scipy.zeros((nN,))
tolerance1 = scipy.zeros((nN,))
tolerance2 = scipy.zeros((nN,))

# creation of three identical liftings with different seeds

sampler_param = inv_transform.Sampler.getDefaultParameters()
sampler_param['seed']=0
lifting0 = inv_transform.Sampler(sampler_param)

sampler_param = inv_transform.Sampler.getDefaultParameters()
sampler_param['seed']=1
lifting1 = inv_transform.Sampler(sampler_param)

sampler_param = inv_transform.Sampler.getDefaultParameters()
sampler_param['seed']=2
lifting2 = inv_transform.Sampler(sampler_param)

for n in range(len(Nlist)):     
    N = Nlist[n]
    param['N']=N
    
    print "running with seed 0 ?"
    fp_sde = particles.Particles(lifting0,restriction,rho,grid,lambd,param)
    rho0_Dt = fp_sde.u_Dt
    tolerance0[n] = scipy.amax(scipy.absolute(rho0_Dt-rho_Dt))
    
    print "running with seed 1 ?" 
    fp_sde = particles.Particles(lifting1,restriction,rho0,grid,lambd,param)
    rho1_Dt = fp_sde.u_Dt
    tolerance1[n] = scipy.amax(scipy.absolute(rho1_Dt-rho_Dt))

    print "running with seed 2 ?" 
    fp_sde = particles.Particles(lifting2,restriction,rho0,grid,lambd,param)
    rho2_Dt = fp_sde.u_Dt
    tolerance2[n] = scipy.amax(scipy.absolute(rho2_Dt-rho_Dt))
