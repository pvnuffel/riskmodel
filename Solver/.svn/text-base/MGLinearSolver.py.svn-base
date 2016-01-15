import LinearSolver
import IdentityLinearSolver

import scipy
from pylab import *

class MGLinearSolver(LinearSolver.LinearSolver):

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
        self.iter = 0
        LinearSolver.LinearSolver.__init__(self,param)
    
    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'smoother':     the name of the smoother
            'presmooth':    the number of presmoothing steps
            'postsmooth':   the number of postsmoothing steps
            'nlevel':       the number of coarsening levels
            'factor':       the coarsening factor 
            'cycle':        the type of multigrid cycle 
                            (at this time, only 'V')
            'nb_iter'       the number of multigrid iterations
            'print'         : 'short'  (one line per iteration)
                other possibilities
                'long' (residual in each iteration)
                'none' (prints nothing)
		'plot' (prints in a format that can be saved to file
			for plotting, i.e. without text explaining)
        """
        param = {}
        param['smoother']   = 'Cheby'
        param['presmooth']  = 2
        param['postsmooth'] = 2
        param['nlevel']     = 4
        param['factor']     = 2
        param['cycle']      = 'V'
        param['nb_iter']    = 10
        param['print']      = 'short'
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def setPoint(self,point):
        self.smoother.setPoint(point)
        LinearSolver.LinearSolver.setPoint(self,point)

    def computeDefect(self,x,rhs,level=0):
        d = rhs - self.point.matvec(x)
        return d

    def restrict(self,d):
        # WATCH OUT : THIS ROUTINE ASSUMES (POWERS OF factor) + 1
        # in the restriction, we want to coarsen the solution
        # vector without throwing away the extra conditions
        # that were added for continuation etc.
        # NOTE THAT COARSENING OF THE EXTRA CONDITIONS
        # IS NOT DEALT WITH YET BY Point.matvec
        # SOLUTION: Oosterlee book, chapter Continuation
        factor = self.param['factor']
        d_grid = d[:len(d)-self.nextra:factor]
        d_c = scipy.r_[d_grid,d[len(d)-self.nextra:]]
        return d_c  

    def prolongate(self,v_c):
        # WATCH OUT : THIS ROUTINE ASSUMES (POWERS OF factor) + 1
        # NOTE THAT COARSENING OF THE EXTRA CONDITIONS
        # IS NOT DEALT WITH YET BY Point.matvec
        # SOLUTION: Oosterlee book, chapter Continuation
        factor = self.param['factor']
        v_c_grid=v_c[:len(v_c)-self.nextra]
        v_grid = scipy.zeros(factor*(len(v_c_grid)-1)+1) 
        f_range = scipy.arange(0.,factor,1.)
        v_grid[::factor]=v_c_grid
        for f in f_range[1:]:
            v_grid[int(f)::factor]=(1.-f/factor)*v_c_grid[:-1]+f/factor*v_c_grid[1:]
        v = scipy.r_[v_grid,v_c[len(v_c)-self.nextra:]]
        return v 
        
    def solve(self,rhs):
        self.nextra = len(rhs)-len(self.point.getState()[0])
        if self.param['cycle']=='V':
            return self.V_solve(rhs)
        elif self.param['cycle']=='W':
            return self.W_solve(rhs)

    def V_solve(self,rhs):
        x = scipy.zeros(len(rhs))
        for i in range(self.param['nb_iter']):
            x, st = self.V_cycle(rhs,x)
            # print statements
            if self.param['print']=='plot':
                res = self.computeDefect(x,rhs)
                res_norm = scipy.linalg.norm(res,2)/len(res)
                print i, res_norm		
	    if self.param['print']=='short'\
                or self.param['print']=='long':
                res = self.computeDefect(x,rhs)
                res_norm = scipy.linalg.norm(res,2)/len(res)
                print "MG: ", i, " | norm res  :", res_norm
            if self.param['print']=='long':
                print "MG Residual"
                for i in range(len(res)):
                    print i, res[i]
                print "---------------------"
        return x, st

    def V_cycle(self,rhs,x,level=0):
        """ 
        Overrides LinearSolver.solve
        Result contains (solution,status)
        Status is always 0; it doesn't mean anything
        x = None :  the current Guess is not given
                    this means it is the current state of the point
            Else:   this is the current guess that needs to be improved
        """
        nlevel = self.param['nlevel']
        presmooth = self.param['presmooth']
        postsmooth = self.param['postsmooth']
        factor = self.param['factor']
       
        if not (level == nlevel):
            # this is when we have not reached the maximum level
            # presmoothing
            for k in range(presmooth):
                x,st = self.smoother.solve(rhs,x)
            # coarse grid correction
            d = self.computeDefect(x,rhs,level)
            d_c = self.restrict(d)
            x_c = scipy.zeros(len(d_c))
            v_c,st = self.V_cycle(factor**2*d_c,x_c,level=level+1)
            v = self.prolongate(v_c)
            x = x + v
            # postsmoothing
            for k in range(postsmooth):
                x,st = self.smoother.solve(rhs,x)
        else:
            # here we have reached the maximum level, so
            # we do a coarse-grid direct solve
            # we are assuming that we only have matvec 
            n = len(rhs)
            A = scipy.zeros((n,n),scipy.float64)
            # we build the (small) system matrix by performing
            # n matvec operations with unit vectors
            for i in range(n):
                vec = scipy.zeros((n,),scipy.float64)
                vec[i]=1.
                A[:,i]= self.point.matvec(vec)
            A = self.point.boundary_conditions(A)
            # and solve using the direct solver of scipy
            x = scipy.linalg.solve(A,rhs) 
        status = 0
        return (x,status)
        
    def W_solve(self,rhs,level=0):
        """ 
        Overrides LinearSolver.solve
        Result contains (solution,status)
        Status is always 0; it doesn't mean anything
        x = None :  the current Guess is not given
                    this means it is the current state of the point
            Else:   this is the current guess that needs to be improved
        """
        if self.param['print']=='long':
            print "Multigrid | level: " , level, " iter: ", self.iter
        nlevel = self.param['nlevel']
        presmooth = self.param['presmooth']
        postsmooth = self.param['postsmooth']
        nb_iter = self.param['nb_iter']

        # current guess
        if level==0:
            x = self.point.getCurrentGuess()
            self.nextra = len(rhs)-len(self.point.getState()[0])
        else:   # we are solving the defect equation on a
                # coarser grid 
            x = scipy.zeros((len(rhs),))
        if not (level == nlevel):
            # this is when we have not reached the maximum level
            for j in range(nb_iter):
                if level == 0:
                    self.iter +=1
                # presmoothing
                for k in range(presmooth):
                    x,st = self.smoother.solve(rhs,x)
                # coarse grid correction
                d = self.computeDefect(x,rhs)
                d_c = self.restrict(d)
                v_c,st = self.W_solve(d_c,level=level+1)
                v = self.prolongate(v_c)
                x = x + v
                # postsmoothing
                for k in range(postsmooth):
                    x,st = self.smoother.solve(rhs,x)
                # print statements
                if level == 0: 
                    if not (self.param['print']=='none'):
                        res = self.computeDefect(x,rhs)
                        res_norm = scipy.linalg.norm(res,2)/len(res)
                        print "MG: ", self.iter, " | norm res  :", res_norm
                    if self.param['print']=='long':
                        print "MG Residual"
                        for i in range(len(res)):
                            print i, res[i]
                        print "---------------------"
        else:
            # here we have reached the maximum level, so
            # we do a coarse-grid direct solve
            # we are assuming that we only have matvec 
            n = len(rhs)
            A = scipy.zeros((n,n),scipy.float64)
            # we build the (small) system matrix by performing
            # n matvec operations with unit vectors
            for i in range(n):
                vec = scipy.zeros((n,),scipy.float64)
                vec[i]=1.
                A[:,i]= self.point.matvec(vec)
            # and solve using the direct solver of scipy
            x = scipy.linalg.solve(A,rhs) 
        status = 0
        return (x,status)
