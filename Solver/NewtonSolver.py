import scipy
import Solver
import LinearSolver
from scipy import zeros,dot,r_

class NewtonSolver(Solver.Solver):

    def __init__(self,linear_solver,parameters=None):
        """ 
        input: 
        =====
            linear_solver (LinearSolver) 
                contains the linear solver that will be used in 
                each Newton iteration
            parameters (dict) 
                look at the docstring of getDefaultParameters() 
                to find out which fields there are
        behaviour:
        =========
            This class implements a Newton solver that stops when the 
            maximum number of iterations has been reached, OR the 
            relative OR absolute tolerance have been reached.
        """
        Solver.Solver.__init__(self,parameters)
        if isinstance(linear_solver,LinearSolver.LinearSolver):
            self.linsolv=linear_solver
        else:
            raise TypeError, "input argument " + linear_solver \
                + " should be a linear solver"
        self.nb_newt = 0 
        self.resid = zeros((0,))
        self.newton_states = zeros((0,))
        self.newton_res_norm = zeros(0,)
       
#    def createCallback(self):
#        def callback(res, x_new):
#            self.resid = r_[self.resid,res]
#          #  self.x = r_[self.x,x_new]
#            if self.param['print']=='short':
#                print "Newton residual: ", len(self.resid), res
#        self.callback=callback
    
    def getDefaultParameters():
        """  
        Returns a dictionary of default parameters 
            'maxIter'       : 10
            'relTol'        : 1e-6
            'absTol'        : 1e-8
            'damping'       : 1 (means no damping)
            'print'         : 'short'  (one line per iteration)
                other possibilities
                'long' (residual in each iteration)
                'none' (prints nothing)
        """
        param = {}
        param['max_iter']=10
        param['rel_tol']=1e-6
        param['abs_tol']=1e-8
        param['print']='short'
        param['damping']=1
        param['stop']="res"
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)
    


    def solve(self):
        """
        Perform Newton iterations, starting from x0.  
        The method returns status
        Status means:
            0: method has converged
            1: method has used maximum number of iterations 
               without convergence
                """
        iter = 0
        max_iter = self.param['max_iter']
        abs_tol = self.param['abs_tol']
        rel_tol = self.param['rel_tol']
        print_it = self.param['print']
        alpha = self.param['damping']
        stop = self.param['stop']
        x = self.point.getCurrentGuess()
        self.newton_states = x
        res = self.point.getResidual()
        res_norm = scipy.linalg.norm(res)/len(res)
        print 'res norm init ', self.newton_res_norm
        self.newton_res_norm  =res_norm
        grid = self.point.system.grid
        grid_dx = grid[-1] - grid[-2]
        if stop == "step":
            stop_norm = 0
        elif stop == "res":
            stop_norm = res_norm
            
        while (iter < max_iter and stop_norm > abs_tol):

            iter += 1
            self.nb_newt += 1
            dx, status = self.linsolv.solve(-res)
            
            step_norm = scipy.linalg.norm(dx)
            step_avg = scipy.sum(dx)
            
            # Handling of errors and printing
            if not status == 0:
                print "An error has occurred in solving " +\
                    "the linear system"
            # print the iterations
            if not print_it=='none':
                print "# Newton: ", iter, " | norm res  :", \
                    res_norm, "norm step : ", step_norm, "avg step: ", step_avg
            if print_it=='long':
                print "Newton Residual"
                for i in range(len(res)):
                    print i, res[i]
                print "---------------------"
                print "Newton Current Guess (before this iteration)"
                for i in range(len(x)):
                    print i, x[i]
                print "---------------------"
                print " Newton Step"
                for i in range(len(x)):
                    print i, dx[i]
                print "---------------------"
                print "Sum : ", scipy.sum(dx)
                print "--------"
                print "Newton Current Guess (after this iteration)"
                for i in range(len(x)):
                    print i, x[i]+alpha*dx[i]
                print "---------------------"
                
            x=x+alpha*dx
          #  x = x/(sum(x)*grid_dx)
            self.newton_states = r_[self.newton_states,x]
            self.point.setCurrentGuess(x)                          
            res = self.point.getResidual()
            res_norm = scipy.linalg.norm(res)/len(res)
            self.newton_res_norm = r_[self.newton_res_norm,res_norm]
            if stop == "step":
                stop_norm = scipy.linalg.norm(dx)
            elif stop == "res":
                stop_norm = res_norm
        if not print_it == 'none':
            print "Newton: final residual : ", res_norm
        #self.iterations = iter
        #self.point.iterations=iter
        if iter == max_iter:
            status = 1
        else:
            status = 0
        return status
        

    def setPoint(self,point):
        self.linsolv.setPoint(point)
        Solver.Solver.setPoint(self,point)
    
    def new_build(self):
        self.linsolv.new_build()       
