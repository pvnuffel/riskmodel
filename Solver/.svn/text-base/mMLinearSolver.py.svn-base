import LinearSolver
import Identity

import scipy

class mMLinearSolver(LinearSolver.LinearSolver):
    """
    This class defines a solver that can be used for
    preconditioning.  It takes a principle from multigrid,
    but the "coarse solve" is not on less grid points; instead
    it is solving the coarse approximate equation for the error
    on the same grid, and lifting the results.
    """

    def __init__(self,param=None):
        """ 
        input: 
        =====
            parameters (dict) 
                You can get a sample dictionary by looking at 
                getDefaultParameters()  
        behaviour:
        =========
            This class implements a Multigrid solver
        """
        if param['smoother']==None:
            self.smoother = Identity.Identity()
        else:
            self.smoother = param['smoother']
        LinearSolver.LinearSolver.__init__(self,param)
    
    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'nlift':        the number of constrained runs steps
                            after the macroscopic solve    
        """
        param = {}
        param['nb_iter'] = 2
        param['additive'] = True
        param['eta'] = 1.
        param['smoother'] = None
        param['presmooth'] = 1
        param['postsmooth'] = 1
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def setPoint(self,point):
        self.smoother.setPoint(point)
        LinearSolver.LinearSolver.setPoint(self,point)

    def computeDefect(self,x,rhs,level=0):
        d = rhs - self.point.matvec(x)
        return d

    def restrict(self,d):
        d_grid = self.point.system.restrict(d[:len(d)-self.nextra])
        d_c = scipy.r_[d_grid,d[len(d)-self.nextra:]]
        return d_c  

    def prolongate(self,v_c):
        # this prolongation does a lifting to micro
        v_c_grid=v_c[:len(v_c)-self.nextra]
        v_grid = self.point.system.lift(v_c_grid)
        v = scipy.r_[v_grid,v_c[len(v_c)-self.nextra:]]
        return v 

    def solve(self,rhs):
        self.nextra = len(rhs)-len(self.point.getState()[0])
        x = scipy.zeros(len(rhs))
        for i in range(self.param['nb_iter']):
            x, st = self.V_cycle(rhs,x)
        return x, st

    def V_cycle(self,rhs,x):
        """ 
        Result contains (solution,status)
        Status is always 0; it doesn't mean anything
        x = None :  the current Guess is not given
                    this means it is the current state of the point
            Else:   this is the current guess that needs to be improved
        """
        
        try:
            for k in range(self.param['presmooth']):
                x,st = self.smoother.solve(rhs,x)
        except: 
            pass
        # take rhs to the coarse grid
        d = self.computeDefect(x,rhs)
        d_c = self.restrict(d)
        x,st = self.smoother.solve(rhs,x)
        # get the coarse Jacobian
        A = self.point.computeCoarseJacobian()
        v_c = scipy.linalg.solve(A,d_c) 
        # back to distributions 
        v = self.prolongate(v_c)
        if self.param['additive']: 
            v = v + self.param['eta'] * d
        x = x + v
        # the other possibility is to lift (x_c + v_c)
        for i in range(self.param['postsmooth']):
            x,st = self.smoother.solve(rhs,x)
        
        status = 0
        
        return (x,status)
