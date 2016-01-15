import scipy
import tables
from scipy.linalg import norm
import TimestepperSystem as TS


class MGTimestepperSystem(TS.TimestepperSystem):

    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and parameter values
        """
        TS.TimestepperSystem.__init__(self,u,lambd,parameters)
        self.computed_n = []
        self.u_Dts = []
   
    def applyJacobian(self,v):
        """
        If the system Jacobian is J, this function returns (an
        approximation to) (1-J)*v 
        """
        eps = self.param['eps']
        u_eps = self.restrict(self.u,len(v)) + eps * v/norm(v)
        u_eps_Dt = self.integrate(u_eps,self.lambd)
        u_Dt = self.get_u_Dt(len(v))
        result = v-(u_eps_Dt - u_Dt)/eps*norm(v)
        return result   

    def restrict(self,u,n):
        factor = int((len(u)-1)/(n-1))
        u_r = u[::factor]
        return u_r

    def get_u_Dt(self,n):
        """ we only want to compute u_Dt on each level once"""
        if n in self.computed_n:
            i = self.computed_n.index(n)
            u_Dt = self.u_Dts[i]
        else:
            self.computed_n.append(n)
            if n == len(self.u):
                u_Dt = self.u_Dt
            else:
                u = self.restrict(self.u,n)
                u_Dt = self.integrate(u,self.lambd)
            self.u_Dts.append(u_Dt)
        return u_Dt
                
