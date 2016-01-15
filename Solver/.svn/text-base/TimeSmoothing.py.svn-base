import LinearSolver

import scipy

class TimeSmoothing(LinearSolver.LinearSolver):

    def __init__(self,parameters=None):
        LinearSolver.LinearSolver.__init__(self,parameters)

    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'tau'       : the Richardson damping factor
            'nb_iter'   : the number of Richardson iterations
            'print'         : 'short'  (one line per iteration)
                other possibilities
                'long' (residual in each iteration)
                'none' (prints nothing)
        """
        param = {}
        param['nb_iter']  = 2
        param['print']  = 'short'
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def solve(self,rhs,x=None,level=0):
        """ 
        Performs a number of Richardson iterations, starting from x
        If no x0 is given, we start from a vector of zeros.
        (If this is used as a stand-alone solver, other routines using
        it, do not expect that an initial guess needs to be given.)
        Returns:
            x :         the Richardson approximation to the solution
            status:     always 0, doesn't mean anything
        """
        status = 0
        if x == None:
            x = scipy.zeros(len(rhs))
        for iter in range(self.param['nb_iter']):
            x = self.point.system.integrate(x,self.point.system.lambd)
        return x,status
