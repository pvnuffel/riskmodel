import types
import tables
import scipy
import scipy.sparse.linalg
from scipy.linalg import norm

import Solver.Solver as Solver
import Solver.LinearSolver as LinearSolver
import System.System as System

class Point(scipy.sparse.linalg.LinearOperator):

    """ 
    This class represents a point structure, which can be corrected.
    As a user, you should only care about the methods
        - getDefaultParameters()
        - getState()/setState()
        - mspectrum(), pspectrum(), spectra()
        - correct()
        - save()

    The following methods are used by the solvers
        - matvec()
        - psolve()
        - getResidual()
        - setCurrentGuess()
        - axpy() (used in continuation)
        
    The following methods are used only internally
        - MatrixVectorProduct()
        - spectrum()
        - getHDF5Description()
    """ 

    def __init__(self,system,solver,precond=None,param=None):
        """ 
        input: 
        =====
            system (System)             
                contains the concrete equations or time-stepper
            solver (Solver) 
                contains the nonlinear solver that will be used 
            precond (LinearSolver) 
                contains the linear solver that will be used as a
                preconditioner
            parameters (dict) 
                look at the docstring of getDefaultParameters() 
                to find out which fields there are
        behaviour:
        =========
            This class implements a Point structure that can be used 
            to store different types of attractors and correct them.
        """
        if isinstance(system,System.System):
            self.system=system
        else:
            raise TypeError, "first argument should be a system"
        if isinstance(solver,Solver.Solver):
            self.solver = solver
            # make sure that the solver knows the point
            self.solver.setPoint(self)
        else:
            raise TypeError, "second argument should be a solver"
        if precond == None:
            self.precond = None
        elif isinstance(precond,LinearSolver.LinearSolver):
            self.precond = precond
            # make sure that the preconditioner knows the point
            self.precond.setPoint(self)
        else:
            raise TypeError, "third argument should be a linear solver"
        if param == None:
            param = self.getDefaultParameters()
        self.param={}
        self.param['free'] = param['free'][:]
        self.param['artificial'] = param['artificial'][:]
        self.param['extra_condition'] = param['extra_condition']
        self.param['artificial_condition'] = param['artificial_condition']
        self.param['sparse_jacobian'] = param['sparse_jacobian'] 
        self.extra_condition = param['extra_condition']
        self.artificial_condition = param['artificial_condition']
        self.nb_matvec=0
        self.u=self.system.u[:]
        self.neq = self.system.neq
        self.lambd=self.system.lambd[:]
        # to be able to use it as a LinearOperator in iterative methods
        self.setShape()
        self.dtype = self.u.dtype
    
    def setShape(self):    
        n = len(self.u)
        l = len(self.param['free'])
        a = len(self.param['artificial'])
        self.shape = (n+l+a,n+l+a)
    
    def getDefaultParameters():
        """  
        Returns a dictionary of default parameters 
            free : list of indices of the free parameters
            artficial  : 
                   list of indices of artificially added free parameters
            extra_condition : a list of functions containing 
                        extra conditions
            artificial_condition : a list of function containing 
                        artificial extra conditions
        """
        param = {}
        param['free']=[]
        param['artificial']=[]
        param['extra_condition']=[]
        param['artificial_condition']=[]
        param['sparse_jacobian']=False
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)
    
    def setState(self,u,lambd):
        self.u=u
        self.lambd=lambd
        self.system.setState(u,lambd)
        # this is meant to deal with preconditioner matrix that needs to be
        # built if the state changes.
        if not self.precond == None:
            self.precond.new_build()
        self.solver.new_build()
    
    def getState(self):
        """
        Returns (u,lambd)
        """
        return self.u,self.lambd
    
    def setCurrentGuess(self,x):
        """
        INPUT:
        =====
        The vector x contains the current guess of the solution.
        
        DATA LAYOUT:
        ===========

        The vector "x" has a part "u"  that describes the
        current solution, "free" for values of the free parameters 
        and a part "art" containing the values of the artificial
        parameters.
        
             [       ]
             [   u   ]
         x = [       ]
             [ ----- ]
             [  free ]
             [  art  ]

        EFFECT:  
        ======
        The current state and solution values are updated.
        """
        n = len(self.u)
        l = len(self.param['free'])
        a = len(self.param['artificial'])
        u = x[:n]
        lambd = self.lambd
        lambd[self.param['free']]=x[n:n+l]
        lambd[self.param['artificial']]=x[n+l:]
        self.setState(u,lambd)
    
    def getCurrentGuess(self):
        """
        OUTPUT:
        =====
        The vector x contains the current guess of the solution.
        
        DATA LAYOUT:
        ===========

        The vector "x" has a part "u"  that describes the
        current solution, "free" for values of the free parameters 
        and a part "art" containing the values of thhe artificial
        parameters.
        
             [       ]
             [   u   ]
         x = [       ]
             [ ----- ]
             [  free ]
             [  art  ]

        """
        u = self.u
        lfree = scipy.take(self.lambd,self.param['free'])
        lart  = scipy.take(self.lambd,self.param['artificial'])
        return scipy.r_[u,lfree,lart]
    
    def getResidual(self):
        """
        Returns the residual associated with the current guess 
        The residual consists of
        (system.getResidual(),residual of extra conditions,
        residual of artificial extra conditions)
        """
        l = len(self.param['free'])
        a = len(self.param['artificial'])
        extra = []
        for i in range(l):
            extra.append(self.extra_condition[i](self))
        artificial = []
        for i in range(a):
            artificial.append(self.artificial_condition[i](self,\
                self.param['artificial'][i]))
        res = self.system.getResidual()
        for i in range(a):
            res +=artificial[i]['eq_term']
        for i in range(l):
            r = extra[i]['res']
            res = scipy.r_[res,r]
        for i in range(a):
            r = artificial[i]['res'] 
            res = scipy.r_[res,r] 
        return res
    
    def computeJacobian(self): 
        """ 
        PURPOSE:
        =======
        This method computes the full system matrix for a linear
        solve.  This consistes of the Jacobian of the underlying
        system, supplemented with the rows and columns that 
        come from additional constraints and free parameters.
        """
        free = self.param['free']
        art = self.param['artificial']
        n=len(self.u)
        l=len(free)
        a=len(art)
        extra = []
        for i in range(l):
            extra.append(self.extra_condition[i](self))
        artificial = []
        for i in range(a):
            artificial.append(self.artificial_condition[i]\
                (self,self.param['artificial'][i]))

        if self.param['sparse_jacobian']:
            J = scipy.sparse.lil_matrix((n+l+a,n+l+a))
        else:
            J = scipy.zeros((n+l+a,n+l+a),scipy.float64)
        A = self.system.computeJacobian()
        if self.param['sparse_jacobian']:
            J[:n,:n] = A.tocsr()
        else:
            J[:n,:n] = A
        # we want the total matrix to be
        #
        #[            |          |  j_art     ]    [   x   ]
        #[  msys      | j_free   | art[column]]    [       ]
        #[            |          |            ]    [       ]
        #[------------------------------------]  * [ ----- ]
        #[ extra[row] | extra[d] |      0     ]    [  free ]
        #[ art[row]   |    0     |    art[d]  ]    [  art  ]
        # 
        J_free = scipy.zeros((n,l),scipy.float64)
        for i in range(l):
            J_free[:,i]=self.system.getParameterDerivative(free[i])
        J[:n,n:n+l]=J_free
        
        J_art = scipy.zeros((n,a),scipy.float64)
        for i in range(a):
            J_art[:,i]=artificial[i]['column']
        
        J[:n,n+l:n+l+a]=J_art
        
        d_free = scipy.zeros((l,l+a),scipy.float64)
        for i in range(l):
            d_free[i,:]=extra[i]['d']
        J[n:n+l,n:n+l+a]=d_free
        
        d_art = scipy.zeros((a,l+a),scipy.float64)
        for i in range(a):
            d_art[i,:]=artificial[i]['d']
        J[n+l:n+l+a,n:n+l+a]=d_art
        
        rows = scipy.zeros((l+a,n),scipy.float64)
        for i in range(l):
            rows[i,:]=extra[i]['row']
        for i in range(a):
            rows[i,:]=artificial[i]['row']
        J[n:n+l+a,:n]=rows

        return J
    
    def boundary_conditions(self,A):
        return self.system.boundary_conditions(A)
    
    def matvec(self,s):
        """
        PURPOSE:
        =======
        Matrix vector product, used in Krylov solvers

        DATA LAYOUT:
        ===========
        The vector "s" has a part "x"  that describes the
        wave and a part "free" for the free parameters and a part for
        the artificial parameters.
        
            [   x   ]
            [       ]
        s=  [ ----- ]
            [  free ]
            [  art  ]
                                

        BACKGROUND:
        ==========
        It uses the Jacobian vector product of the underlying
        system, and adds the effect of extra conditions and free
        parameters.
        The result is a matrix-vector product with the linearization of
        the full nonlinear determining system around the current guess. 
        """
        return self.MatrixVectorProduct(s,"jac")
    
    def pmatvec(self,s):
        """
        PURPOSE:
        =======
        Matrix vector product with the preconditioner matrix, 
        used in Krylov solvers

        DATA LAYOUT:
        ===========
        The vector "s" has a part "x"  that describes the
        wave and a part "free" for the free parameters and a part for
        the artificial parameters.
        
            [   x   ]
            [       ]
        s=  [ ----- ]
            [  free ]
            [  art  ]
                                

        BACKGROUND:
        ==========
        It uses the matrix vector product applyPreconditioner(v) of the underlying
        system, and adds the effect of extra conditions and free
        parameters.
        The result is a matrix-vector product with the linearization of
        the full nonlinear determining system around the current guess. 
        """
        return self.MatrixVectorProduct(s,"precond") 
    
    def MatrixVectorProduct(self,s,option):
        """
        PURPOSE:
        =======
        Matrix vector product, used in Krylov solvers
 
        INPUT:
        =====
        s -- the vector to be multiplied
        option -- 
            "jac"  - compute Jacobian-vector product
            "precond" - preconditioner-vector product

        USE:
        ===
        This method is used in the matvec and pmatvec routines.
        Do not use directly.

        DATA LAYOUT:
        ===========
        The vector "s" has a part "x"  that describes the
        wave and a part "free" for the free parameters and a part for
        the artificial parameters.
        
            [   x   ]
            [       ]
        s=  [ ----- ]
            [  free ]
            [  art  ]
                                

        BACKGROUND:
        ==========
        It uses the matrix vector product applyPreconditioner(v) of the underlying
        system, and adds the effect of extra conditions and free
        parameters.
        The result is a matrix-vector product with the linearization of
        the full nonlinear determining system around the current guess. 
                """
        # setting some local variables
        free=self.param['free']
        art=self.param['artificial']
        l=len(free)
        a=len(art)
        n=len(s)-l-a
        extra = []
        for i in range(l):
            extra.append(self.extra_condition[i](self))
        artificial = []
        for i in range(a):
            artificial.append(self.artificial_condition[i](self,\
                self.param['artificial'][i]))
        sx=s[:n]
        sfree=s[n:n+l]
        sart=s[n+l:]
        sx_norm=norm(sx)
        sfree_norm=norm(sfree)
        
        result=scipy.zeros(scipy.shape(s),scipy.float64)
        
        # obtain the first N rows of the matrix-vector product
        #[        |        |      ]    [   x   ]
        #[ J_x    | J_free |  row ]    [       ]
        #[        |        |      ]    [       ]
        #-------------------------]  * [ ----- ]
        #                              [  free ]
        #                              [  art  ]
        if not sx_norm==0:
            if option == "jac":
                Mv_x=self.system.applyJacobian(sx)
            if option == "precond":
                Mv_x=self.system.applyPreconditioner(sx)
        else:
            Mv_x=sx  # sx will contain zeros
        
        for i in range(l):
            if not sfree[i]==0:
                Mv_x += sfree[i]*self.system.getParameterDerivative(free[i])
        for i in range(a):
            Mv_x += sart[i]*artificial[i]['column']
            
        result[:n]=Mv_x

        sparam = scipy.r_[sfree,sart]
        
        # the rows of the extra conditions are added
        #                              [   x   ]
        #                              [       ]
        #                              [       ]
        #-------------------------]  * [ ----- ]
        #[  row   |   d    |  0   ]    [  free ]
        #                              [  art  ]
        Mv_free=scipy.zeros(scipy.shape(free),scipy.float64)
        for i in range(l):
            Mv_free[i]=scipy.dot(extra[i]['row'],sx)+\
                    scipy.dot(extra[i]['d'],sparam)
        result[n:n+l]=Mv_free
        # the rows of the artificial conditions are added
        #                              [   x   ]
        #                              [       ]
        #                              [       ]
        #-------------------------]  * [ ----- ]
        #                              [  free ]
        # [ row    | 0      | d   ]    [  art  ]
        Mv_art=scipy.zeros(scipy.shape(art),scipy.float64)
        for i in range(a):
            Mv_art[i]=scipy.dot(artificial[i]['row'],sx)+\
                    scipy.dot(artificial[i]['d'],sparam)
        result[n+l:]=Mv_art
        return result
    
    def correct(self):
        """
        Corrects the initial guess using self.solver.
        Returns 0 upon convergence
        """
        return self.solver.solve()
    
    def psolve(self,rhs):
        """ 
        PURPOSE:
        =======
        This method solves a linear system M*x = b, which 
        closely ressemles the linear system A*x = b as a
        preconditioner.
        This solve is done together with each matrix-vector product
        inside any Krylov method.
        The matrix M contains the preconitioner matrix of the system,
        bordered with the extra and artificial conditions.
        
        INPUT:
        =====
        rhs     Right-hand side
        
        """
        if not self.precond == None:
            x,status = self.precond.solve(rhs)
        else:
            x = rhs.copy()
        return x
    
    def mspectrum(self):
        """ 
        Returns the spectrum of the Jacobian
        OUTPUT:
        ======
            e -- eigenvalues of the Jacobian
            v -- eigenvectors of the Jacobian
        """
        return self.spectrum("jac")
    
    def pspectrum(self):
        """ 
        Returns the spectrum of the preconditioner
        OUTPUT:
        ======
            e -- eigenvalues of the Jacobian
            v -- eigenvectors of the Jacobian
        """
        return self.spectrum("precond")
    
    def spectrum(self,option):
        """ 
        Returns the spectrum of the Jacobian or the preconditioner.
        Note: this is used internally.
        
        INPUT:
        ======
            option -- 
                "jac"  - compute Jacobian spectrum
                "precond" - compute preconditioner spectrum

        OUTPUT:
        ======
            e -- eigenvalues of the Jacobian
            v -- eigenvectors of the Jacobian
        """
        n=len(self.u)
        l=len(self.param['free'])
        a=len(self.param['art'])
        A=scipy.zeros((n+l+a,n+l+a),scipy.float64)
        for i in range(n+l+a):
            y=scipy.zeros((n+l+a,),scipy.float64)
            y[i]=1
            # matrix of original system
            if option == "jac":
                A[i,:]=self.matvec(y)
            if option == "precond":
                A[i,:]=self.psolve(y)
        eA,vA = scipy.linalg.eig(A)
        if option == "precond":
            eA = 1./eA
        return eA,vA
    
    def spectra(self):
        """ 
        Returns the spectrum of the Jacobian A, the preconditioner M,
        and the resulting matrix M^(-1)*A .
        
        OUTPUT:
        ======
            ej -- eigenvalues of the Jacobian
            ep -- eigenvalues of the preconditioner
            ec -- eigenvalues of the product matrix
            vj -- eigenvectors of the Jacobian
            vp -- eigenvectors of the preconditioner
            vc -- eigenvectors of the product matrix
        """
        n=len(self.u)
        l=len(self.param['free'])
        a=len(self.param['artificial'])
        A=scipy.zeros((n+l+a,n+l+a),scipy.float64)
        B=scipy.zeros((n+l+a,n+l+a),scipy.float64)
        C=scipy.zeros((n+l+a,n+l+a),scipy.float64)
        for i in range(n+l+a):
            y=scipy.zeros((n+l+a,),scipy.float64)
            y[i]=1
            # matrix of original system
            A[i,:]=self.matvec(y)
            # matrix of preconditioned system
            B[i,:]=self.psolve(A[i,:])
            # matrix of preconditioner -- inverse of it
            C[i,:]=self.psolve(y)
        ej,vj = scipy.linalg.eig(A)
        ec,vc = scipy.linalg.eig(B)
        ep,vp = scipy.linalg.eig(C)
        ep = 1./ep
        return ej,ep,ec,vj,vp,vc        
    
    def axpy(self,a,pointy):
        """
        Returns a * point + pointy
        """
        p=Point(self.system,self.solver,self.param)
        p.setState(a*self.u+pointy.u,\
            a*self.lambd+pointy.lambd)
        return p
    
    def copy(self):
        return Point(self.system,self.solver,self.param)
    
    def save(self,file,group,table):
        # om zeker te zijn dat het systeem in een
        # coherente state gesaved wordt!
        self.system.setState(self.u,self.lambd)
        self.system.save(file,group,table)
    
    def getHDF5description(self):
        return self.system.getHDF5Description()
    
