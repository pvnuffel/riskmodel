from scipy.sparse.linalg import gmres as sgmres
#from iterative import gmres as sgmres
#from iterative import bicgstab as bicgstab

from scipy.linalg import norm, lstsq
from scipy import zeros,dot,r_
import LinearSolver

class GMRESLinearSolver(LinearSolver.LinearSolver):

    def __init__(self,param=None):
        """ 
        input: 
        =====
            parameters (dict) 
                You can get a sample dictionary by looking at 
                getDefaultParameters()  
        behaviour:
        =========
            This class implements a GMRES solver that stops when the 
            maximum number of iterations has been reached, OR the 
            relative OR absolute tolerance have been reached.
        """
        LinearSolver.LinearSolver.__init__(self,param)
        self.resid = zeros((0,))
        self.createCallback()
    
    def createCallback(self):
        def callback(res):
            self.resid = r_[self.resid,res]
            if self.param['print']=='short':
                print "gmres residual: ", len(self.resid), res
        self.callback=callback
    
    def getDefaultParameters():
        """
        Returns a dictionary of default parameters 
            'print'      : 'none' (prints nothing)
                other possible values
                    'short' (print a line for every iteration)
                    'long' (prints solution)
            'restart'    : number of iterations before GMRES restarts
            'tol'        : relative tolerance to achieve
            'relative'   : False
                If True, the 'tol' parameter is considered relative
                with respect to the right-hand side.
        """
        param = {}
        param['print']='short'
        param['tol']= 1e-5
        param['relative'] = False
        param['rel_tol']=0.5
        param['max_iter']= 5000
        param['restart']= None
        param['flexible']= False
        param['builtin']= True
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def solve(self,rhs):
        """ 
        Overrides LinearSolver.solve
        Result contains (solution,status)
        Status means:
             0: method has converged
             k>0 : maximum number of iterations reached
             k<0 : illegal input or breakdown
        """
        self.resid = zeros((0,))
        tol = self.param['tol']
        if self.param['relative']:
            tol = max(self.param['rel_tol']*norm(rhs),tol)
        print_it = self.param['print']
        #if print_it == 'short':
        restart = self.param['restart']
        maxiter = self.param['max_iter']
        if self.param['flexible']:
            x,info=self.fgmres(rhs,tol=tol,restrt=restart,\
                maxiter=maxiter)
        elif self.param['builtin']:
            if ((restart > len(rhs)) or (restart == None)):
                restart = 2*len(rhs)
            print "# using builtin | restrt : ", restart
            x,info=sgmres(self.point,rhs,tol=tol,\
                restrt=restart,maxiter=maxiter,\
                callback=self.callback,xtype=None)
        else:
            x,info=self.gmres(rhs,tol=tol,restrt=restart,\
                maxiter=maxiter)
        return (x,info)


    def fgmres(self,rhs,tol=1e-6,restrt=None,maxiter=None,callback=None):
        if maxiter == None:
            maxiter = len(rhs)
        if restrt == None:
            restrt = 2*maxiter
        # implemented as in [Saad, 1993]
        # start
        x = zeros(len(rhs))
        H = zeros((restrt+1, restrt))
        V = zeros((len(rhs),restrt))
        Z = zeros((len(rhs),restrt))
        # Arnoldi process (with modified Gramm-Schmidt)
        res = 1.
        j = 0
        r = rhs - self.point.matvec(x)
        beta = norm(r)
        V[:,0]=r/beta
        while j < maxiter and res > tol:
            Z[:,j] = self.point.psolve(V[:,j])
            w = self.point.matvec(Z[:,j])
            for i in range(j+1):
                H[i,j]=dot(w,V[:,i])
                w = w - H[i,j]*V[:,i]
            H[j+1,j] = norm(w)
            V[:,j+1]=w/H[j+1,j]
            e = zeros(j+2)
            e[0]=1.
            y, res, rank, sing_val = lstsq(H[:j+2,:j+1],beta*e)
            j += 1
            print "# GMRES| iteration :", j, "res: ", res/beta
            self.resid = r_[self.resid,res/beta]
            Zy = dot(Z[:,:j],y)
        x = x + Zy
        info = 1
        return (x,info)

    def gmres(self,rhs,tol=1e-6,restrt=None,maxiter=None,callback=None):
        if maxiter == None:
            maxiter = len(rhs)
        if restrt == None:
            restrt =  maxiter
        # implemented as in [Saad, 1993]
        # start
        n = len(rhs)
        x = zeros(n)
        H = zeros((n+1, n))
        V = zeros((n,n))
        Z = zeros((n,n))
        # Arnoldi process (with modified Gramm-Schmidt)
        j = 0
        r = rhs - self.point.matvec(x)
        beta = norm(r)
        myres = beta/len(r)
        V[:,0]=r/beta
        print "in gmres",
        print myres,tol
        while j < maxiter and myres > tol:
            z = self.point.psolve(V[:,j])
            w = self.point.matvec(z)
            for i in range(j+1):
                H[i,j]=dot(w,V[:,i])
                w = w - H[i,j]*V[:,i]

            H[j+1,j] = norm(w)
            V[:,j+1]=w/H[j+1,j]
            e = zeros(j+2)
            e[0]=1.
            y, res, rank, sing_val = lstsq(H[:j+2,:j+1],beta*e)
            print "in my own gmres: ", j, myres
            self.resid = r_[self.resid,myres]
            j += 1
            Vy = dot(V[:,:j],y)
            xbar = x + self.point.psolve(Vy)
            rbar = rhs - self.point.matvec(xbar)
            myres = norm(rbar)/len(rbar)
        try: 
            x = x + self.point.psolve(Vy)
        except UnboundLocalError:
            pass
        info = 0
        return (x,info)
