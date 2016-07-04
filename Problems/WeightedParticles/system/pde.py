import scipy
import numpy as np


import sys         
sys.path.append('../')   #added to resolve importerror of nephew-directories
sys.path.append('../../')
sys.path.append('../../../')

import System.TimestepperSystem as TSystem

from scipy.integrate import quad

import matplotlib.pylab as plt

def doublewell(x):
    return 2*x-4*x**3   

def d2udx2(u,dx):
    u_bc=scipy.r_[u[0],u,u[-1]]        #add first and last element
    #(so they appear two times and u_bc[1:-1] is original vector with same dimension then u_bc[:-2] )
    return 1./dx**2*(u_bc[2:]-2*u_bc[1:-1]+u_bc[:-2]) #FTCS scheme

    
class FokkerPlanck(TSystem.TimestepperSystem):
    """
    Implements a PDE based time-stepper for the PDE
    
    rho_t + [a*(x-x**3) \rho]_x = D \rho_xx 
    """
    def __init__(self,rho,grid,drift,lambd,param=None):
        self.grid = grid
        xL = 3./2.*grid[0]-1./2.*grid[1]
        xR = 3./2.*grid[-1]-1./2.*grid[-2]
        self.domain = scipy.array([xL,xR])
        self.drift = drift
        if param == None :
            param = FokkerPlanck.getDefaultParameters()
        self.param = param
        self.neq = 1
      #  self.precond = param['precond']
 #       self.precond.system = self
        self.Dt = param['Dt']
        self.dt = param['dt'] 
        TSystem.TimestepperSystem.__init__(self,rho,lambd,param)    
    
    def integrate(self,rho0,lambd):
        Dt = self.param['Dt']
        # Dt needs to be a multiple of param['Dt']
        dt = self.param['dt']
        tcur = 0.
        rho = scipy.array(rho0)
        plot_rho_dt = True
        n_steps = int(Dt/dt)
        for tcur in range(0, n_steps):
            rho = rho + dt*self.rhs(rho,lambd)
            if(plot_rho_dt and tcur % 100 == 0):
                self.makePlot(rho,lambd, tcur)
        return rho
    
    def rhs(self,rho,lambd):
        grid = self.grid
        dx = grid[1]-grid[0]
        D = lambd[1]
        a = lambd[0]
        flux = self.flux(rho)
        rhs = -a*(flux[1:]-flux[:-1])/dx+D*d2udx2(rho,dx)
        return rhs

    def flux(self,rho):
        edges = scipy.r_[self.domain[0],(self.grid[:-1]+self.grid[1:])/2.,self.domain[-1]]
        flux = scipy.zeros_like(edges)
        flux[1:-1]=self.drift(edges[1:-1])*(rho[1:]+rho[:-1])/2.
        return flux
    
    def computeJacobian(self,precond=False):
        if precond:
            dx = self.grid[1]-self.grid[0]
            N = len(self.u)
            return self.precond.computeJacobian(N,dx) 
        else:
            return TSystem.TimestepperSystem.computeJacobian(self)
    
    def makePlot(self, rho, lambd, tcur):   
        plt.axis([-2, 2, 0, 0.5]) 
        label1= r'$\rho (x,t)$'
        plot_rho = plt.plot(self.grid, rho, label=label1)
        plt.ylabel(r'$\rho$', fontsize = 16)
        plt.xlabel('$x$', fontsize = 14)
        D = lambd[1]
        a =lambd[0]
        dx = self.grid[1]-self.grid[0]
        norm = sum( np.exp((-self.grid**4 + self.grid**2 )*a/D))*dx
   #         norm = quad(lambda x: np.exp(-( x**4 -x**2)*a/D ), -np.inf, np.inf)  #gives same_result
        rho_ss = np.exp( (-self.grid**4 + self.grid**2 )*a/D)/norm
        label2= r'$\frac{ \exp{\left[-\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}$'
        plot_rho = plt.plot(self.grid, rho_ss, label=label2)
        plt.legend(prop={'size':0.1})
        plt.legend([plot_rho], loc='best')
        plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
        plt.savefig('movieplots/plot_rho_t%.5d.jpg' %tcur, dpi=500)
        #plt.savefig('movieplots/plot_rho_t%.5d.pdf' %tcur)
        plt.close()      
    
