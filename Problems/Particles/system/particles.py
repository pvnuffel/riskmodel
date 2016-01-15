import scipy

import sys         
sys.path.append('../')   #added to resolve importerror of nephew-directories
sys.path.append('../../')
sys.path.append('../../../')




import System.TimestepperSystem as TSystem

def doublewell(x):
    return x-x**3

class Particles(TSystem.TimestepperSystem):
    """
    Implements a particle-based time-stepper for an ensemble of particles
    evolving according to 
    dx = a dt + sqrt(2*D)*dW
    
    The corresponding Fokker-Planck is 
    
    rho_t + [a*(x-x**3) \rho]_x = D \rho_xx 
    """
    def __init__(self,lifting,restriction,rho,grid,lambd,param=None):
        print param
        if param == None :
            param = Particles.getDefaultParameters()
        self.param = param
        self.grid = grid
        # grid points are cell centers
        xL = 3./2.*grid[0]-1./2.*grid[1]
        xR = 3./2.*grid[-1]-1./2.*grid[-2]
        self.domain = scipy.array([xL,xR])
        self.lifting = lifting
        self.restriction = restriction
        self.neq = 1
        self.rand=scipy.random.mtrand.RandomState()
        #self.precond = param['precond']
        #self.precond.system = self
        self.Dt = param['Dt']
        TSystem.TimestepperSystem.__init__(self,rho,lambd,param)
    
    def getDefaultParameters():
        param = {}
        param['drift']=doublewell
        param['N']=10000
        param['dt']=1e-3
        param['Dt']=1e-2
        param['seed']=0
        param['eps']=1e-5
        return param 
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def seed(self,s):
        self.cur_seed = s
        self.rand.seed(self.cur_seed)
    
    def integrate(self,rho,lambd):
        # lambd = [a,D]
        self.lifting.seed(self.lifting.param['seed'])
        #self.seed(self.param['seed'])
        x = self.lift(rho)
        print "nb particles : ", len(x)
        x_Dt = self.simulate(x,lambd)
        return self.restrict(x_Dt)
    
    def lift(self,rho):
        grid = self.grid
        x = self.lifting.lift(rho,grid,self.param['N'])
        return x
    
    def simulate(self,x0,lambd):
        Dt = self.param['Dt']
        # Dt needs to be a multiple of param['Dt']
        dt = self.param['dt']
        D = lambd[1]
        a = lambd[0]
        N = self.param['N']
        drift = self.param['drift']
        x = scipy.array(x0)
        
        tstart = 0
        tcur = tstart
        while (tcur < tstart + Dt + dt/2 ):
            tcur += dt
            # the random number
            dW=self.rand.normal(loc=0.,scale=scipy.sqrt(2*D*dt),size=N)
          #  if tcur == dt:    #only print random number for first time step
          #      print 'dW =',  dW            
            
            # the process
            drift_term = a * drift(x)
            x=x+drift_term*dt+dW
            # and reflecting boundary conditions
            scipy.where(x>self.domain[1],2*self.domain[1]-x,x)
            scipy.where(x<self.domain[0],2*self.domain[0]-x,x)
        return x
    
    def restrict(self,x):
        return self.restriction.restrict(x,self.grid,self.domain)
    
    def computeJacobian(self):
        dx = self.grid[1]-self.grid[0]
        N = len(self.u)
        return self.precond.computeJacobian(N,dx)
    
