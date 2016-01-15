import scipy
import tables
from scipy.linalg import norm
import TimestepperSystem as TS


class mMTimestepperSystem(TS.TimestepperSystem):

    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and parameter values
        """
        TS.TimestepperSystem.__init__(self,u,lambd,parameters)
        
    def setState(self,u,lambd):
        TS.TimestepperSystem.setState(self,u,lambd)
        self.U = self.restrict(u)
        self.U_Dt = self.restrict(self.u_Dt)

    def restrict(self,u):
        """
        A restriction operator
        """
        assert 1, "Needs to be implemented for each system"
        
    
    def lift(self,U):
        """
        A restriction operator
        """
        assert 1, "Needs to be implemented for each system"