import tables
import scipy
import sys

from scipy.linalg import norm

import matplotlib.pylab as plt

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

import inv_transform_sampling as inv_transform
 #import inv_transform_determ as inv_transform  (same results for different seeds)
import histogram
import particles
import pde
import precond
#import TimestepperSystem

if __name__=="__main__":
     
    D = None
    Dt = None
    N = None
    M=1  #number of monte carlo steps
    # deze manier van input parameters verwerken staat in Langtangen
    while len(sys.argv) > 1 :
        option = sys.argv[1];               del sys.argv[1]
        if option == "-Dt" :
            Dt = float(sys.argv[1]);         del sys.argv[1]
            print "#Dt: ", Dt, " (from command line)"
        elif option == "-N" :
            N = int(sys.argv[1]);         del sys.argv[1]
            print "#N: ", N, " (from command line)"
            
    if D == None:
        D = 1./2.
    if Dt == None:
        Dt = 1e-2
        
    if N == None:
        Nlarge = 1000
        Nsmall = 1000
    seed = 1
    # discretization parameters
    h=2e-2 # kde mesh width 
    dx = 0.01
    dt = 1e-5 # ti
    r= (2.*D*dt)/dx**2
    if r>1: 
        print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
        sys.exit("Ending Script. Make sure that r<1") 
    xL = -1.7
    xR = 1.7
    grid = scipy.arange(xL+dx/2.,xR,dx)

    rho = scipy.ones_like(grid)
    #rho[150]=1
        # 0.1*scipy.random.uniform(low=-0.5,high=0.5,size=scipy.shape(grid)[0])
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho

    # Fokker Planck for SDE
    a = 1
    zeta = 0.
  #  alphas = [1.]
    alpha=1
    beta = -1

    p = None
    fp_sde = None
    

    
   # for i in range(len(alphas)):
    
  #      alpha = alphas[i]
#        print "alpha ... ",  alpha

        # Fokker Planck for SDE
    lambd = scipy.array([a,D,alpha,beta,zeta])

    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=seed
    lifting = inv_transform.Sampler(sampler_param)
    
    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=0
    lifting0 =  inv_transform.Sampler(sampler_param)

    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=1
    lifting1 = inv_transform.Sampler(sampler_param)
   
    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=h
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
    param['Dt'] = Dt
    param['Nlarge']=Nlarge
    param['Nsmall']=Nsmall   
    param['dt']=dt
    precond_param=precond.Precond.getDefaultParameters()
    precond_param['Dstar']=1./2.
    precond_param['sigma']=scipy.zeros_like(grid)
    precond_param['kappa']=scipy.zeros_like(grid)
    param['precond']=precond.Precond(precond_param)

    param_histogram = histogram.Histogram.getDefaultParameters()
    restriction = histogram.Histogram(param_histogram)
    
    fprev = fp_sde
    fp_sde = particles.Particles(lifting,restriction,rho,grid,\
        lambd,control = fprev, param=param)
    print fp_sde,fp_sde.control
    # print fp_sde2,fp_sde2.control
    
    fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
 
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
    newt_param['max_iter']=1
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

    #v = scipy.ones_like(p.u)
    v=scipy.zeros_like(p.u)
    v[20]=-1
    v[40]=1

    # u_eps_Dt,u_eps_Dt_control,result=\
    #     fp_sde.testJacobian(v)  
    
    Jv_sde=fp_sde.applyJacobian(v)
    rho_Dt = fp_sde.u_Dt
    Jv_pde = fp_pde.applyJacobian(v)    
    rho_Dt_fp = fp_pde.u_Dt
    norm_Jv_fp = norm(Jv_pde)
    

        
#     
#    print "running with seed 0 ?"
#    fp_sde=None
#    fp_sde = particles.Particles(lifting0,restriction,rho,grid,lambd,fp_sde,param)
#    rho0_Dt = fp_sde.u_Dt
#    tolerance0 = scipy.amax(scipy.absolute(rho0_Dt-rho_Dt_fp))
#    Jv_sde0 =fp_sde.applyJacobian(v)     
#    
#    
#    print "running with seed 1 ?" 
#    fp_sde=None
#    fp_sde = particles.Particles(lifting1,restriction,rho,grid,lambd,fp_sde,param)
#    rho1_Dt = fp_sde.u_Dt
#    tolerance1 = scipy.amax(scipy.absolute(rho1_Dt-rho_Dt_fp))
#    Jv_sde1 =fp_sde.applyJacobian(v)  
#
#    print "running with seed 2 ?" 
#    fp_sde=None
#    fp_sde = particles.Particles(lifting2,restreiction,rho,grid,lambd,fp_sde,param)
#    rho2_Dt = fp_sde.u_Dt
#    tolerance2 = scipy.amax(scipy.absolute(rho2_Dt-rho_Dt_fp))
#    Jv_sde2 =fp_sde.applyJacobian(v) 
    
    Nlist = scipy.array([1000,2000,5000,10000])     
    
    tolerance = scipy.zeros((M+1))
    Jv_norm = scipy.zeros((M+1))
    diff_norm = scipy.zeros((M+1))
    rho_norm = scipy.zeros((M+1))
    rho_42 = scipy.zeros((M+1))
    Jv_sde_sum = 0
    totaldiff = 0
    rho42=0
    rho_average=0
    #rho_m =scipy.zeros_like(rho_Dt)
    
    
    
    for m in range(1,1):
        print "running with seed ", m 
        sampler_param = inv_transform.Sampler.getDefaultParameters()
        sampler_param['seed']=m
        lifting = inv_transform.Sampler(sampler_param)            
                
        fp_sde=None
        fp_sde = particles.Particles(lifting,restriction,rho,grid,lambd,fp_sde,param)
        rho_Dt = fp_sde.u_Dt
#        rho42=(rho42+ rho_Dt[42])/m
#        rho_42[m] = rho42 -  rho_Dt_fp[42]
#        rho_average = (rho_average + rho_Dt)
#        rho_norm[m] = norm((rho_average - rho_Dt_fp)/m)      \
        
     #   tolerance[m] = scipy.amax(scipy.absolute(rho_Dt-rho_Dt_fp))
        Jv_sde =fp_sde.applyJacobian(v)          
        Jv_sde_sum = Jv_sde_sum + Jv_sde
        diff =  Jv_sde - Jv_pde
        totaldiff = totaldiff + diff
        averagediff = totaldiff/(m)
        diff_norm[m] = norm(averagediff)
        Jv_norm[m] = norm(Jv_sde_sum/(m))
        
        #norm(Jv_sde-Jv_pde)  #twonorm of vector difference
       
      
       # plt.plot(u_eps_Dt)
  