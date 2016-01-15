import LinearSolver

import scipy

class Cheby(LinearSolver.LinearSolver):
    """
    T_k(x) is minimal norm in [-1,1].
    We want minimal norm in [l_max/2,l_max] with constraint T_k(0)=1
    We therefore take alpha*T_k((x-3*l_max/4)/(l_max/2)). Alpha is the scaling
    to ensure T(0)=1 
    """
    def __init__(self,parameters=None):
        LinearSolver.LinearSolver.__init__(self,parameters)

    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'lambd_max'     : the maximal eigenvalue for the problem 
                we do optimal smoothing for [lambd_max/2..lambd_max] 
            'order'         : the order of the polynomial
                default = 2
            'nb_iter'   : the number of iterations
            'print'         : 'short'  (one line per iteration)
                other possibilities
                'long' (residual in each iteration)
                'none' (prints nothing)
        """
        param = {}
        param['lambd_max'] = 1.
        param['order'] = 2
        param['nb_iter']  = 1
        param['print']  = 'short'
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def polyvec(self,res):
        lmax = self.param['lambd_max']
        order = self.param['order']
        coeff = scipy.zeros(order)
        if order == 2:
            coeff[1]=-32./lmax**2/17.
            coeff[0]=48./lmax/17.
        else:
            print "for this order, the coefficients have not been coded"
        t_1 = self.point.matvec(res)
        return coeff[1]*t_1+coeff[0]*res

    def solve(self,rhs,x=None):
        """ 
        Performs a number of polynomial smoothing iterations, 
        starting from x
        If no x0 is given, we start from a vector of zeros.
        (If this is used as a stand-alone solver, other routines using
        it do not expect that an initial guess needs to be given.)
        Returns:
            x :         the Chebychev approximation to the solution
            status:     always 0, doesn't mean anything
        """
        status = 0
        if x == None:
            x = scipy.zeros(len(rhs))
        for iter in range(self.param['nb_iter']):
            res = rhs - self.point.matvec(x)
            x = x + self.polyvec(res)
            if not self.param['print']=='none':
                res_norm = scipy.linalg.norm(res,2)/len(res) 
                print "Cheby: ", iter, " | norm res  :", res_norm
            if self.param['print']=='long':
                print "Cheby Residual"
                for i in range(len(res)):
                    print i, res[i]
                print "---------------------"
        return x,status
