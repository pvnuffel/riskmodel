import scipy
import scipy.sparse.linalg
import LinearSolver
        
class ImplicitEulerDirectLinearSolver(LinearSolver.LinearSolver):
    """
    This linear solver is written, based on the derivation of the implicit
    Euler simplification as indicated in G. Samaey (PhD, p. 127 bottom).
    It assumes that the Jacobian that is returned from the linear
    system is given by
    Dt * J -- and not (I - Jacobian of timestepper)
    As a consequence, the rhs and the added rows and columns need to be
    adjusted accordingly.
    """
    
    def __init__(self,param=None):
        LinearSolver.LinearSolver.__init__(self,param)
        self.built=False
    
    def new_build(self):
        self.built=False
        
    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'print'      : 'none' (prints nothing)
                other possible value
                    'long' (prints solution)
        """
        param = {}
        param['print']='none'
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def solve(self,rhs):
        """ 
        Overrides LinearSolver.solve
        Result contains (solution,status)
            status is always 0, indicating that the method has converged
        """
        rhs_ = rhs.copy()
        N = len(self.point.getState()[0])
        Dt = self.point.system.Dt
        I = scipy.identity(N,scipy.float64)
        
        if not self.built:
            A = self.point.computeJacobian()
            cols = A[:,N:]
            if type(A)==scipy.sparse.lil_matrix:
                cols-=Dt*A*cols
                A[:,N:]=-cols
            else:
                c = scipy.shape(cols)[1]
                for i in range(c):
                    cols[:N,i]=-scipy.dot(I-Dt*A[:N,:N],cols[:N,i])
                A[:,N:] = cols
            self.A = A
            self.built = True
        else:
            A = self.A
                    
        if type(A)==scipy.sparse.lil_matrix:
            rhs_ = -rhs_+Dt*A*rhs_
            rhs_[-1] = rhs[-1]
            x=scipy.sparse.linalg.spsolve(Dt*A.tocsr(),rhs_)
        else:
            rhs_[:N] = -scipy.dot(I-Dt*A[:N,:N],rhs_[:N])
            x=scipy.linalg.solve(Dt*A,rhs_)
        
        status = 0
        return (x,status)
        
