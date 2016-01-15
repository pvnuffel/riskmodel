import scipy

import EquationSystem as EqS

class  MGEquationSystem(EqS.EquationSystem):

    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and 
        parameter values
        """
        EqS.EquationSystem.__init__(self,u,lambd,parameters)
    
    def computeJacobian(self,n=None):
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
        A = self.computeJacobian(len(v))
        return scipy.dot(A,v)
