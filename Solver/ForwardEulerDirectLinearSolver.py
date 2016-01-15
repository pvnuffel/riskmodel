import scipy
import LinearSolver
        
class ForwardEulerDirectLinearSolver(LinearSolver.LinearSolver):
    
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
        if not self.built:
            N = len(self.point.getState()[0])
            Dt = self.point.system.Dt
            dt = self.point.system.dt
            k = int(Dt/dt)
            I = scipy.identity(N)
            
            A = self.point.computeJacobian()
            B = I+dt*A[:N,:N]
            AA = B
            for i in range(k-1):
                AA=scipy.dot(AA,B)
            Matrix = A  # zo blijven extra rijen en kolommen dezelfde als A
            Matrix[:N,:N]= I - AA
            self.Matrix = Matrix
            self.built = True
        else:
            Matrix=self.Matrix
        x=scipy.linalg.solve(Matrix,rhs)
        status = 0
        return (x,status)
        
