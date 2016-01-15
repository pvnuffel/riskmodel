import LinearSolver

import scipy

class Richardson(LinearSolver.LinearSolver):

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
        param['tau']   = 1.
        param['nb_iter']  = 2
        param['print']  = 'short'
        param['tol'] = 1e-8
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
        tau = self.param['tau']*(level+1)**2
        max_iter = self.param['nb_iter']
        tol = self.param['tol']
        status = 0
        if x == None:
            x = scipy.zeros(len(rhs))
        res_norm = 1.
        it = 0
        while ((it < max_iter) and (res_norm > tol)):
            res = rhs - self.point.matvec(x)
            prhs = self.point.psolve(rhs)
            # x = x + tau * self.point.psolve(res)
            x = x - tau * self.point.psolve(self.point.matvec(x))\
                + tau * prhs                    
            if not self.param['print']=='none':
                res_norm = scipy.linalg.norm(res,2)/len(res) 
                print it, res_norm, scipy.linalg.norm(x,2)/len(x)
            if self.param['print']=='long':
                print "Richardson Residual"
                for i in range(len(res)):
                    print i, res[i]
                print "---------------------"
            it += 1
        return x,status
