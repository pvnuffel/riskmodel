import LinearSolver

import scipy

class Identity(LinearSolver.LinearSolver):

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
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def solve(self,rhs,x=None):
        """ 
        This solver just returns the rhs.
        Its only use is to be able to simulate unpreconditioned GMRES easily.
        """
        status = 0
        x = scipy.array(rhs)
        return x,status
