import types
import tables
import scipy
from scipy.linalg.iterative import gmres
from scipy.linalg import norm

import Point
import Solver.Solver as Solver
import Solver.LinearSolver as LinearSolver
import System.System as System

class mMPoint(Point.Point):

    """ 
    This class represents a mMpoint structure, which can be corrected.
    It inherits from Point
    """ 

    def __init__(self,system,solver,precond=None,param=None):
        """ 
        input: 
        =====
            system (System)             
                contains the concrete equations or time-stepper
            solver (Solver) 
                contains the nonlinear solver that will be used 
            precond (LinearSolver) 
                contains the linear solver that will be used as a
                preconditioner
            parameters (dict) 
                look at the docstring of getDefaultParameters() 
                to find out which fields there are
        behaviour:
        =========
            This class implements a Point structure that can be used 
            to store different types of attractors and correct them.
        """
        Point.Point.__init__(self,system,solver,precond,param)
    
    def getDefaultParameters():
        """  
        Returns a dictionary of default parameters 
            free : list of indices of the free parameters
            artficial  : 
                   list of indices of artificially added free parameters
            extra_condition : a list of functions containing 
                        extra conditions
            artificial_condition : a list of function containing 
                        artificial extra conditions
        """
        param = {}
        param['free']=[]
        param['artificial']=[]
        param['extra_condition']=[]
        param['artificial_condition']=[]
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def computeCoarseJacobian(self): 
        """ 
        PURPOSE:
        =======
        This method computes the full system matrix for a linear
        solve.  This consistes of the Jacobian of the underlying
        system, supplemented with the rows and columns that 
        come from additional constraints and free parameters.
        """
        # since this is an mMpoint, we know that the system has a self.u and a
        # self.U; here we want to compute with the self.U state, but everything
        # else is expecting the state to be in self.u
        self.u= self.system.U
        self.neq = self.system.Neq
    
        free = self.param['free']
        art = self.param['artificial']
        n=len(self.u)
        l=len(free)
        a=len(art)
        extra = []
        for i in range(l):
            extra.append(self.extra_condition[i](self))
        artificial = []
        for i in range(a):
            artificial.append(self.artificial_condition[i]\
                (self,self.param['artificial'][i]))

        J = scipy.zeros((n+l+a,n+l+a),scipy.float64)
        J[:n,:n] = self.system.computeCoarseJacobian()

        # we want the total matrix to be
        #
        #[            |          |  j_art     ]    [   x   ]
        #[  msys      | j_free   | art[column]]    [       ]
        #[            |          |            ]    [       ]
        #[------------------------------------]  * [ ----- ]
        #[ extra[row] | extra[d] |      0     ]    [  free ]
        #[ art[row]   |    0     |    art[d]  ]    [  art  ]
        # 
        J_free = scipy.zeros((n,l),scipy.float64)
        for i in range(l):
            J_free[:,i]=self.system.getParameterDerivative(free[i])
        J[:n,n:n+l]=J_free
    
        J_art = scipy.zeros((n,a),scipy.float64)
        for i in range(a):
            J_art[:,i]=artificial[i]['column']
        
        J[:n,n+l:n+l+a]=J_art
        
        d_free = scipy.zeros((l,l+a),scipy.float64)
        for i in range(l):
            d_free[i,:]=extra[i]['d']
        J[n:n+l,n:n+l+a]=d_free
        
        d_art = scipy.zeros((a,l+a),scipy.float64)
        for i in range(a):
            d_art[i,:]=artificial[i]['d']
        J[n+l:n+l+a,n:n+l+a]=d_art
        
        rows = scipy.zeros((l+a,n),scipy.float64)
        for i in range(l):
            rows[i,:]=extra[i]['row']
        for i in range(a):
            rows[i,:]=artificial[i]['row']
        J[n:n+l+a,:n]=rows
        
        self.u = self.system.u
        self.neq = self.system.neq

        return J

    