import scipy
from scipy.linalg import norm
from scipy import sqrt
import numpy as np
from numpy import *
import sys         
sys.path.append('../')   #added to resolve importerror of nephew-directories
sys.path.append('../../')
sys.path.append('../../../')

import System.TimestepperSystem as TSystem

def create_doublewell(alpha,beta):
    def doublewell(x):
       # return -4*x**3+3*(alpha+beta)*x**2-2*alpha*beta*x
      #  return 2*x-4*x**3  
        return x-x**3  
    return doublewell

class Particles(TSystem.TimestepperSystem):
    """
    Implements a particle-based time-stepper for an ensemble of particles
    evolving according to 
    dx = a dt + sqrt(2*D)*dW
    
    The corresponding Fokker-Planck is 
    
    rho_t + [a*(x-x**3) \rho]_x = D \rho_xx 
    """
    def __init__(self,lifting,restriction,rho,grid,lambd, x_prev=None, w_prev=None, control=None,param=None):
   
        if param == None :
            param = Particles.getDefaultParameters()
        self.param = param
        self.drift=create_doublewell(alpha=lambd[2],beta=lambd[3])
        self.control= control
#        if control == None:
#            self.N=param['Nlarge']
#            self.skip=1
#        else:
#            self.N=param['Nsmall']
#            self.skip = param['Nlarge']/param['Nsmall']
        self.skip=1
        self.N= param['Nsmall']
        self.Nlarge= param['Nlarge']
        self.grid = grid
        # grid points are cell centers
        xL = 3./2.*grid[0]-1./2.*grid[1]
        xR = 3./2.*grid[-1]-1./2.*grid[-2]
        self.domain = scipy.array([xL,xR])
       # self.discr_factor= param['discr_factor']
        self.lifting = lifting
        self.restriction = restriction
        self.neq = 1
        self.rand=scipy.random.mtrand.RandomState() #https://docs.scipy.org/doc/numpy-1.6.0/reference/generated/numpy.random.mtrand.RandomState.html
     #   if(rg_state != None):
   #         self.rand.set_state(rg_state) #tried this for new random numbers in simulations with different particles numbers, but not needed

  #      self.precond = param['precond']    #  (  keyError: 'precond' )
#       self.precond.system = self
        self.Dt = param['Dt']
        self.dt = param['dt']
        self.x_prev = x_prev
        self.w_prev = w_prev
       # print "begin constructor Tsystem"
        TSystem.TimestepperSystem.__init__(self,rho,lambd,param)
        #  print "end constructor"
    
    def getDefaultParameters():
        param = {}
        param['Nlarge']=1000000
        param['Nsmall']=10000
        param['dt']=1e-4
        param['Dt']=1e-2
        param['seed']=0
        param['eps']=1e-5
        return param 
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def seed_particle(self,s):
        self.rand.seed(s) 
        #If seed is None, then RandomState will try to read data from /dev/urandom 
        #if available or seed from the clock otherwise.
    
    def setState(self,u,lambd):
        self.u = u
        self.lambd = lambd
        self.x,self.x_Dt,self.u_Dt, self.x_mean = self.integrate(u,lambd)
        self.setBins()
    
    def integrate(self,rho,lambd):
        # lambd = [a,D]
     #   print "begin lifting"
        x= self.lift(rho)
     #   print "end lifting"
 #       print "x= ", x
       # print "nb particles : ", len(x)
        #seed = self.rand.seed(self.param['seed']) 
      #  seed =self.lifting.seed(self.lifting.param['seed'])  #=number of realization
        seed =self.lifting.param['seed']
        seedlist = scipy.zeros_like(x)      #We don't use this seeds anymore!
        for i in range(len(seedlist)):
            seedlist[i] = seed*self.Nlarge + i  #fill list with  same elements (need this datatype for inv_transform_sampling.py)
     #   print seedlist
       # print "seed= ", seed
        print "begin simulation"
        x_Dt, x_mean = self.simulate(x, lambd)
        if(self.x_prev != None):
            x_Dt = np.concatenate([self.x_prev, x_Dt])
            #print ' x_prev = ', self.x_prev, ' concatenated to ', x_Dt 
        #self.x_prev = x_Dt #this update is already done in loop over particles 
       # else: x_Dt_tot = x_Dt 
        print "restrict"
        rho_new= self.restrict(x_Dt)         #geen gewichten
        return x,x_Dt, rho_new, x_mean  #x_Dt might have other dimension than x
    
    def lift(self,rho):
        grid = self.grid
        x = self.lifting.lift(rho,grid,self.param['Nsmall'],self.skip)
        return x
    
#    def simulate(self,x0,seedlist,lambd):
#        #print "in simulate : ", lambd
#        x_Dt = scipy.zeros_like(x0)
#        for i in range(len(x0)):
#            x_Dt[i]=self.simulate_particle(x0[i],seedlist[i],lambd, mean(x0))
##            if (i==42): 
##                print 'position diff for particle 42 =' , x_Dt[i]-x0[i]
#        return x_Dt
    
#    def simulate_particle(self,x,seed,lambd, mean):
#        tstart = 0
#        tcur = tstart
#        Dt = self.param['Dt']
#        # Dt needs to be a multiple of param['Dt']
#        dt = self.param['dt']
#        D = lambd[1]
#        a = lambd[0]
#        alpha =lambd[2]
#       # print "seed", seed
#        print 'SIMULATE'
#      #  self.seed_particle(seed)    
#        while (tcur < tstart + Dt - dt/2 ):
#            tcur += dt
#            # the random number
#            dW = self.getBrownianIncrement()
#            if tcur == dt:    #only print random number for first time step
#                print dW
#            # the process
#            drift_term = a * self.drift(x)
#            x=x+ (drift_term + alpha*(mean-x))*dt+sqrt(2*D*dt)*dW 
#            # and reflecting boundary conditions
#            if (x>self.domain[1]):
#                x = 2*self.domain[1]-x
#            if (x<self.domain[0]):
#                x = 2*self.domain[0]-x
#        return x
#        
        
        
    def simulate(self,x, lambd):
        Dt = self.param['Dt']
        dt = self.param['dt']
        sigma = lambd[1]
        a = lambd[0]
        alpha =lambd[2]
        
        n_steps = int(Dt/dt)
        x_mean = scipy.zeros(n_steps)
        
        for tcur in range(0, n_steps):
            # the random number
            dW=self.rand.normal(loc=0.,scale=sigma*scipy.sqrt(dt),size=self.N)
           # if tcur == 1:    #only print random number for first time step
            #    print 'dW =',  dW            
                
                    # the process
            drift_term = a * self.drift(x)
            x=x+(drift_term+alpha*(mean(x)-x)) *dt+dW 
            # and reflecting boundary conditions
            scipy.where(x>self.domain[1],2*self.domain[1]-x,x)
            scipy.where(x<self.domain[0],2*self.domain[0]-x,x)
            x_mean[tcur]= mean(x) 
        return x, x_mean
#    
#    
#    
    def getBrownianIncrement(self):
        dW=self.rand.normal(loc=0.,scale=1.)
        return dW
    
    def restrict(self,x):
        return self.restriction.restrict(x,self.grid,self.domain)
        
    def restrict_on_grid(self,x, newgrid):     #niet nodig, want grid lijkt enkel belangrijk voor restrictie
        xL = 3./2.*newgrid[0]-1./2.*newgrid[1]
        xR = 3./2.*newgrid[-1]-1./2.*newgrid[-2]
        newdomain = scipy.array([xL,xR])        
        return self.restriction.restrict(x, newgrid, newdomain)
        
    def make_rho_fine(self,rho):             #obsolete (in stead of making the pde solution fine, we make the pde solution coarse now)
        rho_fine = scipy.zeros(len(rho)*self.discr_factor)
        for i in range (0,len(rho)):
            for j in range (0, self.discr_factor):
                rho_fine[i*self.discr_factor+j]=rho[i]
        return rho_fine                
                        
    
    def setBins(self):   #convert position list for every particle in a list of bin-numbers for ervery particle
        self.bin = self.restriction.getBins(self.x,self.grid,self.domain)
    
    def computeJacobian(self,precond=False):
        if precond:
            dx = self.grid[1]-self.grid[0]
            N = len(self.u)
            return self.precond.computeJacobian(N,dx) 
        else:
            return TSystem.TimestepperSystem.computeJacobian(self)
    
#    def applyJacobian_obs(self,v):
#        eps = self.param['eps']
#        # compute weights per bin for perturbation 
#        # w = 1 + eps * v/rho
#        rho_fine = self.make_rho_fine(self.u) 
#        w_bin = 1. + eps * v/(norm(v)*rho_fine)
#        # transform to weights per particle
#        w_part = w_bin[self.bin]
#        total = scipy.sum(w_part)/self.N
#        # note that the perturbed state is no longer a probability distribution
#        # so norming the histograms is not entirely correct
#        u_eps_Dt = total*self.restriction.restrict(self.x_Dt,self.grid,self.domain,w=w_part)
#        rho_fine_eps_Dt = self.make_rho_fine(u_eps_Dt) 
#        rho_fine_Dt = self.make_rho_fine(self.u_Dt) 
#        Jv = v-(rho_fine_eps_Dt - rho_fine_Dt)/eps*norm(v)
#        #control = self.controlJacobian(v)
#        return Jv  #, u_eps_Dt, self.u_Dt
#        
        
    def applyJacobian(self,v):
        eps = self.param['eps']
        # compute weights per bin for perturbation 
        # w = 1 + eps * v/rhoT
        w_bin = 1. + eps * v/(norm(v)*self.u)
        self.wbin = w_bin
        if(self.w_prev != None):
            bin_tot =  scipy.r_[self.w_prev, self.bin]  
        else:
            bin_tot= self.bin
     #   self.w_prev= bin_tot  #UNCOMMENTED BECAUSE OF P.CORRECT() - STILL NEED TO HANDLE THIS PROBLEM
      #  print 'bintot =',  bin_tot  
        # transform to weights per particle
        w_part = w_bin[bin_tot]   #gives the weight corresponding to the bin index for each particle
   #     print 'Len w_part', len(w_part)
        total = scipy.sum(w_part)/self.Nlarge  #may not be exactly 1
       # print total
        # note that the perturbed state is no longer a probability distribution
        # so norming the histograms is not entirely correct
        self.x_Dt
    #    print 'Len x', len(self.x_Dt)
        u_eps_Dt = total*self.restriction.restrict(self.x_Dt,self.grid,self.domain,w=w_part) 
        #the difference between U_eps_Dt and U_Dt is w_part (they use the same x_Dt)
      #  Jv = v-(u_eps_Dt -self.u_Dt)/eps*norm(v)
        Jv = (u_eps_Dt -self.u_Dt)/eps*norm(v)
        Jv = v-Jv
       # test_same_random_numbers = self.restriction.restrict(self.x_Dt, self.grid, self.domain) -self.u_Dt #=0?
        #control = self.controlJacobian(v)
        return Jv #, u_eps_Dt, self.u_Dt
    
    
    def testJacobian(self,v):
        eps = self.param['eps']
        # compute weights per bin for perturbation 
        # w = 1 + eps * v/rho
        w_bin = 1. + eps * v/norm(v)/self.u
        # transform to weights per particle
        w_part = w_bin[self.bin]
        total = scipy.sum(w_part)/self.N   
        # note that the perturbed state is no longer a probability distribution
        # so norming the histograms is not entirely correct
        u_eps_Dt = total*self.restriction.restrict(self.x_Dt,self.grid,self.domain,w=w_part)
        result = v-(u_eps_Dt - self.u_Dt)/eps*norm(v)
        control = self.controlJacobian(v)
        return result + control,result,control     
    
    def controlJacobian(self,v):
        eps = self.param['eps']
        # compute weights per bin for perturbation 
        # w = 1 + eps * v/rho
        w_bin = 1. + eps * v/norm(v)/self.u
        # transform to weights per particle
        w_part = w_bin[self.bin]
        total = scipy.sum(w_part)/self.N #divide by Nlarge
        # note that the perturbed state is no longer a probability distribution
        # so norming the histograms is not entirely correct
        u_eps_Dt = total*self.restriction.restrict(self.x_Dt,self.grid,self.domain,w=w_part)
        c = self.control
        # print "self.control", c,
        # print self.lambd
        if c==None:
            print "no control variable ..."
            control = scipy.zeros_like(v)
        else: 
            print "control variable ..."
            eps = self.param['eps']
            skip = self.skip
            Nlarge = self.param['Nlarge']
            Nsmall = self.param['Nsmall']
            w_part_large = scipy.ones_like(c.x)
            w_part_small = scipy.ones_like(self.x)
            w_bin_eps = 1. + eps * v/norm(v)/c.u
            w_eps_part_large = w_bin[c.bin]
            w_eps_part_small = w_bin[c.bin[skip/2:Nlarge:skip]]
            total_large = scipy.sum(w_part_large)/Nlarge
            total_small = scipy.sum(w_part_small)/Nsmall
            u_Dt_large = total_large*self.restriction.restrict(\
                c.x_Dt,c.grid,c.domain,w=w_part_large)
            u_eps_Dt_large = total_large*self.restriction.restrict(\
                c.x_Dt,c.grid,c.domain,w=w_eps_part_large)
            control_large = u_eps_Dt_large - u_Dt_large
            u_Dt_small = total_small*self.restriction.restrict(\
                c.x_Dt[skip/2:Nlarge:skip],c.grid,c.domain,w=w_part_small)
            u_eps_Dt_small = total_small*self.restriction.restrict(\
                c.x_Dt[skip/2:Nlarge:skip],c.grid,c.domain,w=w_eps_part_small)
            control_small = (u_eps_Dt_small - u_Dt_small)
            control = (control_small - control_large)/eps*norm(v)
        return control      

    
    
