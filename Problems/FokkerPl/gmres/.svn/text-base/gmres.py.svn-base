import scipy
import numpy
from math import sqrt
import string
import sys
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../..")

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ForwardEulerDirectLinearSolver as FEDirect
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import system.system as system
import system.log as syslog
import system.precond  as precond
import Utils.conditions.periodic as periodic

def beta(x):
    return x-x**3

def D(x):
    return 0.126
    
if __name__=="__main__":
    # p = None
    # 
    # # deze manier van input parameters verwerken staat in Langtangen
    # while len(sys.argv) > 1 :
    #     option = sys.argv[1];               del sys.argv[1]
    #     if option == "-p" :
    #         p = float(sys.argv[1]);         del sys.argv[1]
    # 
    # if p == None:
    #     p = 1.
    
    L = 4.
    dx = 0.05
    x = scipy.arange(-L/2.+dx/2.,L/2,dx)
    y = scipy.exp(-x**2)
    y = y/sum(y)
    u = scipy.array(y)
    
    lambd = scipy.array([0.])
   
    dt = 1e-3
    nsteps = 10
    # CREATING THE SYSTEM OBJECT
    # Changing the default parameters
    param = system.System.getDefaultParameters()
    param['eps'] = 1e-7 
    param['dt'] = dt
    param['x'] = x
    nsteps = 1
    param['Dt'] = nsteps * dt
    param['D']=D
    param['beta']=beta
    precond_param=None
    param['precond']=precond.Precond(precond_param)

    syst = system.System(u,lambd,param)
  
    # CREATING LINEAR SOLVER
    gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
    gmres_param['tol']=2e-6
    gmres_param['print']='short'
    gmres_param['max_iter']=100
    gmres_param['builtin']= True
    linsolv = GMRES.GMRESLinearSolver(gmres_param)

    # CREATING NEWTON SOLVER
    newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
    newt_param['rel_tol']=1e-10
    newt_param['abs_tol']=1e-10
    newt_param['print']='short'
    newt_param['max_iter']=5
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)

    # CREATING POINT
    point_param = Point.Point.getDefaultParameters()
    # point_param['artificial_condition']=[periodic.condition]
    # point_param['artificial']=scipy.array([0])
    psolver = ImDirect.ImplicitEulerDirectLinearSolver()

    print "computing the solution with the linear equations"
    p = Point.Point(syst,nsolv,psolver,point_param)
    p.correct()
    
    print "---------------------"
    
    print "computing the solutions with the log equations"
    # ulog = scipy.log(p.u)+0.1
    ulog = scipy.log(u)
    print "log-solution: ", scipy.log(p.u)
    print "log-initial : ", ulog
    print "difference  : ", scipy.log(p.u)-ulog, scipy.linalg.norm(scipy.log(p.u)-ulog)
    syst_log = syslog.System(ulog,lambd,param)
    plog = Point.Point(syst_log,nsolv,psolver,point_param)
    plog.correct()
    for i in range(len(p.u)):
        print x[i],p.u[i],scipy.exp(plog.u[i])
