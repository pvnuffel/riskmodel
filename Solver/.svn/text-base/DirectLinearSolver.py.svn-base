import scipy
from scipy.sparse import lil_matrix as lil_matrix
import scipy.sparse.linalg.dsolve as dsolve

import LinearSolver
        
class DirectLinearSolver(LinearSolver.LinearSolver):
    
    def __init__(self,param=None):
        LinearSolver.LinearSolver.__init__(self,param)
        self.built=False       
 
    def new_build(self):
        self.built=False       
 
    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'print'      : 'none' (prints nothing)
                other possible value
                    'long' (prints solution)
        """
        param = {}
        param['print']='none'
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def solve(self,rhs):
        """ 
        Overrides LinearSolver.solve
        Result contains (solution,status)
            status is always 0, indicating that the method has converged
        """
        if not self.built:
            A = self.point.computeJacobian()
            self.A=A
        else:
            A=self.A
        if type(A)==lil_matrix:
            x = dsolve(A,rhs)
        else:
            x=scipy.linalg.solve(A,rhs)
        status = 0
        return (x,status)
        
