import tables
import scipy
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

import inv_transform
import histogram
import particles
import pde
import precond

if __name__=="__main__":
     
    D = None
    Dt = None
    N = None
    # deze manier van input parameters verwerken staat in Langtangen
    while len(sys.argv) > 1 :
        option = sys.argv[1];               del sys.argv[1]
        if option == "-D" :
            D = float(sys.argv[1]);         del sys.argv[1]
            print "#D: ", D, " (from command line)"
        elif option == "-Dt" :
            Dt = float(sys.argv[1]);         del sys.argv[1]
            print "#Dt: ", Dt, " (from command line)"
        elif option == "-N" :
            N = int(sys.argv[1]);         del sys.argv[1]
            print "#N: ", N, " (from command line)"
    
    if D == None:
        D = 0.01
    if Dt == None:
        Dt = 0.01 
    if N == None:
        N = 100000
     
    # discretization parameters
    h=2e-2 # kde mesh width 
    dx = 0.05
    dt = 1e-3 
    print "nu : ", dx**2/2./D/dt
    xL = -2.
    xR = 2.
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)/(xR-xL)
    print "rho: ", rho

    # Fokker Planck for SDE
    a = 1
    D = 0.1
    lambd = scipy.array([a,D])

    Dt = 0.01

    # Fokker Planck for SDE
    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=0
    lifting = inv_transform.Sampler(sampler_param)
   
    param_histogram = histogram.Histogram.getDefaultParameters()
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
    param['Dt'] = Dt
    param['N']=100000    
    precond_param=precond.Precond.getDefaultParameters()
    precond_param['Dstar']=0.
    precond_param['sigma']=scipy.zeros_like(grid)
    precond_param['kappa']=particles.doublewell
    param['precond']=precond.Precond(precond_param)

    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=2e-2
    restriction = histogram.Histogram(param_kde)

    fp_sde = particles.Particles(lifting,restriction,rho,grid,lambd,param)

    # CREATING LINEAR SOLVER
    gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
    gmres_param['tol']=1e-8
    gmres_param['print']='short'
    gmres_param['builtin']=True
    linsolv = GMRES.GMRESLinearSolver(gmres_param)

    # CREATING NEWTON SOLVER
    newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
    newt_param['rel_tol']=1e-7
    newt_param['abs_tol']=1e-7
    newt_param['print']='none'
    newt_param['max_iter']=1
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)

    # CREATING POINT
    
    psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()

    p = Point.Point(fp_sde,nsolv,psolver_im)
    p.correct()
        
 
