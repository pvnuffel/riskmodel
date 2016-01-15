import scipy

class Precond(object):
    def __init__(self,param):
        self.Dstar = param['Dstar']
        self.kappa = param['kappa']
        self.sigma = param['sigma']
        self.system = None
    
    def getDefaultParameters():
        param = {}
        param['sigma']=0.
        param['Dstar']=1.
        param['kappa']=0.
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)
    
    def computeJacobian(self,N,dx):
        system = self.system
        Dstar = self.Dstar
        kappa = self.kappa
        sigma = self.sigma 

        edges = (self.system.grid[:-1]+self.system.grid[1:])/2.
        a=kappa(edges)
                
        J=scipy.zeros((N,N),scipy.float64)

        for i in range(1,N-1):
            J[i,i-1] =Dstar/dx**2 + a[i-1]/(2*dx)      
            J[i,i]   =-2*Dstar/dx**2 + sigma[i]-(a[i]-a[i-1])/(2*dx)
            J[i,i+1] =Dstar/dx**2 - a[i]/(2*dx)     
                                  
        # Homogeneous Dirichlet boundary conditions with boundary points eliminated
        
        # left boundary
        J[0,0]   =-2*Dstar/dx**2 + sigma[0] -a[0]/(2*dx)
        J[0,1] =Dstar/dx**2 - a[0]/(2*dx)
        
        # right boundary
        J[N-1,N-2] = Dstar/dx**2 + a[N-2]/(2*dx)  
        J[N-1,N-1]   = -2*Dstar/dx**2 + sigma[N-1] + a[N-2]/(2*dx)
        return J
