import scipy

import System

class EquationSystem(System.System):
    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and 
        parameter values
        """
        System.System.__init__(self,u,lambd,parameters)
    
    def getDefaultParameters():
        return {}
    getDefaultParameters = staticmethod(getDefaultParameters)
        
    def getResidual(self):
        """
        This method determines what the residual for the system is.
        For a fixed point, the method should return a vector containing
        the right-hand side of the system equations
        """
        raise NotImplementedError, \
            "EquationSystem.getResidual() " +\
            "needs to be implemented for a specific system"
        
    def computeJacobian(self):
        """
        returns the Jacobian of the system for the given state and
        parameter vectors
        """
        raise NotImplementedError, \
            "EquationSystem.computeJacobian() " +\
            "needs to be implemented for a specific system"

    def applyJacobian(self,v):
        """
        Return the Jacobian-vector product of the
        system Jacobian with the vector v
        """
        A = self.computeJacobian()
        return scipy.dot(A,v)

    def solvePreconditioner(self,rhs,rows,cols):
        """ 
        solves a linear system with the preconditioning matrix
        
        input:
        =====
            rhs     contains the right-hand side of the system to solve
            rows    contains a number of extra rows that come from external
                constraints
            cols    contains a number of extra columns that contain entries 
                stemming from free parameters
        """
        raise NotImplementedError, \
            "EquationSystem.solvePreconditioner() " +\
            "needs to be implemented for a specific system"
    
    def applyPreconditioner(self,v):
        """ 
        Return a matrix vector product with the preconditioning matrix
        """
        raise NotImplementedError, \
            "EquationSystem.applyPreconditioner() " +\
            "needs to be implemented for a specific system"
        
    def getParameterDerivative(self,i):
        """
        This method returns the derivative with respect to parameter i
        """
        raise NotImplementedError, \
            "EquationSystem.getParameterDerivative() " +\
            "needs to be implemented for a specific system"
