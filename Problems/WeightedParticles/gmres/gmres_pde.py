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
import Utils.conditions.probability as probability

import inv_transform
import histogram
import particles
import pde
import precond

if __name__=="__main__":
     
    D = None
    Dt = None
    # deze manier van input parameters verwerken staat in Langtangen
    while len(sys.argv) > 1 :
        option = sys.argv[1];               del sys.argv[1]
        if option == "-D" :
            D = float(sys.argv[1]);         del sys.argv[1]
            print "#D: ", D, " (from command line)"
        elif option == "-Dt" :
            Dt = float(sys.argv[1]);         del sys.argv[1]
            print "#Dt: ", Dt, " (from command line)"
    
    if D == None:
        D = 1./2.
    if Dt == None:
        Dt = 0.01 
     
    # discretization parameters
    dx = 0.05
    dt = 1e-3 
    print "nu : ", dx**2/2./D/dt
    xL = -1.7
    xR = 1.7
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)/(xR-xL)
    print "rho: ", rho

    # Fokker Planck for PDE
    a = 1
    alpha = 0.99
    lambd = scipy.array([a,D,alpha])

    param=pde.FokkerPlanck.getDefaultParameters()
    param['Dt'] = Dt
    precond_param=precond.Precond.getDefaultParameters()
    precond_param['Dstar']=1./2.
    precond_param['sigma']=scipy.zeros_like(grid)
    # precond_param['kappa']=pde.doublewell
    precond_param['kappa']=scipy.zeros_like(grid)
    param['precond']=precond.Precond(precond_param)

    fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
    print fp_pde.param,fp_pde.Dt
    
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
    newt_param['print']='short'
    newt_param['max_iter']=3
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)

    # CREATING POINT
    
    psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()

    # POINT PARAMETERS
    point_param = Point.Point.getDefaultParameters()
    point_param['artificial']=[2]
    point_param['artificial_condition']=[probability.condition]

    p = Point.Point(fp_pde,nsolv,None,point_param)
    print "before correct: ", sum(p.u)
    p.correct()
    print "after correct: ", sum(p.u)
    
     