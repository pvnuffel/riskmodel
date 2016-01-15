import scipy

class Precond(object):
    def __init__(self,lambd):
        self.lambd = lambd
        self.system = None
    
    def computeJacobian(self):
        system = self.system
        dx = system.dx
        Nx = len(self.system.x)
        
        mu =  1./dx**2
        sigma = -10
        
        J=scipy.zeros((Nx,Nx),scipy.float64)
        for k in range(1,Nx-1):
            # J[k,k-1] = mu
            # J[k,k+1] = mu
            # J[k,k]   =-2*mu+sigma
            J[k,k]=1.
            
        # Homogeneous Dirichlet boundary conditions with boundary points eliminated
        # J[0,1]   = mu
        # J[0,0]   = -2*mu + sigma
        # J[0,-1] = mu
        # J[-1,-2] = mu
        # J[-1,-1] = -2*mu + sigma
        # J[-1,0] = mu
        J[0,0]=1.
        J[-1,-1]=1.
        return J
    

        