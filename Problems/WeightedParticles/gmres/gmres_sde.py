import tables
import scipy
import numpy
import sys

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../..")

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability

from scipy.linalg import norm
from scipy import sqrt

#import inv_transform_sampling as inv_transform
#import inv_transform_determ as inv_transform
import inv_transform_sampling as inv_transform
import histogram
import particles
import pde
import precond


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
    alpha= 1
    beta = -1
    lambd = scipy.array([a,sigma,alpha,beta,zeta])
            
    #SDE
    t1 = time.time()
    
    print 'Start solving sde'
    dx = 1e-2
#    factor = int (dx/dx_fp)
#    print "Discretisation for solving sde is ", factor, " times coarser than the discretisation for solving the pde"
    dt = 1e-2 # ti
    N=1000
    Nlarge= N
    Nsmall =N
    
    
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
   # rho_left = scipy.ones(len(grid)/2)
   # rho_right= scipy.ones(len(grid)/2)/100
   # rho =scipy.r_[rho_left, rho_right] 
    rho = rho/(sum(rho)*dx)
    print rho

    
    h=1e-2
    
    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))

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
    param['seed']=3
    
    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=1
    lifting = inv_transform.Sampler(sampler_param)   
    

#    precond_param=precond.Precond.getDefaultParameters()
#    precond_param['Dstar']=1./2.
#    precond_param['sigma']=scipy.zeros_like(grid)
#    precond_param['kappa']=scipy.zeros_like(grid)
#    param['precond']=precond.Precond(precond_param)
    
    fprev = fp_sde
    x_prev_sim=None
    w_prev_sim=None
    print '!'
    fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, x_prev_sim, w_prev_sim, param=param )

    #_sde = particles.Particles(lifting,restriction,rho,grid,   lambd,control = fprev, param=param)
    print fp_sde
    # print fp_sde2,fp_sde2.control

    # CREATING LINEAR SOLVER
    gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
    gmres_param['tol']=1e-8
    gmres_param['print']='short'
    gmres_param['builtin']=True
    linsolv = GMRES.GMRESLinearSolver(gmres_param)
    # linsolv2 = GMRES.GMRESLinearSolver(gmres_param)

    # CREATING NEWTON SOLVER
    newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
    newt_param['rel_tol']=1e-7
    newt_param['abs_tol']=1e-7
    newt_param['print']='short'
    newt_param['max_iter']=4
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
    # nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)

    # CREATING POINT
    psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
    # psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolver()

    # POINT PARAMETERS
    point_param = Point.Point.getDefaultParameters()
    point_param['artificial']=[4]
    point_param['artificial_condition']=[probability.condition]
    
    pprev = p
    p = Point.Point(fp_sde,nsolv,psolver_im,point_param)
    # p2 = Point.Point(fp_sde2,nsolv2,psolver_im2,point_param)

    p.correct()
