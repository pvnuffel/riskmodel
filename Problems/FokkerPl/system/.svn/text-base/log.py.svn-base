import scipy
import tables
import string
from scipy.linalg import norm

import System.TimestepperSystem as TSystem

class System(TSystem.TimestepperSystem):
    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and parameter values
        Input parameters:
        ================
            state : current system state as a vector
            lambd : current parameters as a vector [D c artificial ?? alpha]
            mesh  : the mesh on which the state is defined
            parameters : see getDefaultParameters() for a list
        """
        self.neq = 1
        if parameters == None :
            parameters = self.getDefaultParameters()
        self.dt = parameters['dt']
        self.x = parameters['x']
        self.dx = self.x[1]-self.x[0]
        self.Dt = parameters['Dt']
        self.D  = parameters['D']
        self.beta = parameters['beta']
        self.precond = parameters['precond']
        self.precond.system = self
        TSystem.TimestepperSystem.__init__(self,u,lambd,parameters)
    
    def getResidual(self):
        """
        This method determines what the residual for the system is.
        """
        res = scipy.exp(self.u) - scipy.exp(self.u_Dt)
        # res = self.u - self.u_Dt
        print "residual : ", norm(res), scipy.sum(res)
        print "res: ", res
        print "+++"
        return res

    def applyJacobian(self,v):
        """
        If the system Jacobian is J, this function returns (an
        approximation to) (1-J)*v 
        """
        result = TSystem.TimestepperSystem.applyJacobian(self,v) 
        result = result * scipy.exp(self.u)   
        result = result/norm(result)
        return result
        
    def rhs(self,u,lambd):
        """
        Dit is de vergelijking
        u_t =  (D(x) u)_xx  +  (beta(x) * u)_x
        """
        dx = self.dx
        flux = self.flux(u)
        rhs = -(flux[1:]-flux[:-1])/dx 
        return rhs
    
    def flux(self,u):
        x = self.x
        dx = self.dx
        Du = self.D(x)*u
        n = scipy.shape(x)[0]
        flux = scipy.zeros((n+1,))
        bu = self.beta(x)*u
        diff = (Du[1:]-Du[:-1])/dx
        adv =scipy.where(self.beta(x[:-1]+dx/2.)>0,bu[:-1],bu[1:])
        flux[1:-1]= adv - diff
        return flux
    
    def integrate(self,u,lambd):
        dt = self.dt
        Dt = self.Dt
        dx = self.dx
        alpha, steps = scipy.modf(Dt/dt)
        # we assume that u is the logarithm of the density
        # print "in integrate"
        # print "u: ", u
        # print "+++"
        y=scipy.exp(u)
        # print "y: ", y
        # print "+++"
        for i in range(int(steps)):
               rhs=self.rhs(y,lambd)
               y=y+rhs*dt
        rhs=self.rhs(y,lambd)
        y_new=y+rhs*dt
        y_Dt = alpha*y_new+(1-alpha)*y
        # print "y_Dt: ", y_Dt
        # print '++++'
        # print sum(y_Dt)
        # print "+++"
        u_Dt = scipy.log(y_Dt)
        # print "u_Dt: ", u_Dt
        # print "-------"
        return u_Dt
    
    def computeJacobian(self):
        return self.precond.computeJacobian() 
    
    def getHDF5Description(self):
        u_shape = scipy.shape(self.u)
        mesh_shape = scipy.shape(self.mesh)
        lambd_shape = scipy.shape(self.lambd)
        class TSystemSave(tables.IsDescription):
            u = tables.FloatCol(shape = u_shape)
            mesh = tables.FloatCol(shape = mesh_shape)
            lambd = tables.FloatCol(shape = lambd_shape)
            dt = tables.FloatCol()
            Dt = tables.FloatCol()
            eps = tables.FloatCol()
            dx = tables.FloatCol()
        return TSystemSave
        
    
    def save(self,file,group,table):
        h5file = tables.openFile(file,mode="a")
        table = self.getTable(h5file,group,table)
        point=table.row
        point['u']=self.u
        point['lambd']=self.lambd
        point['dt']=self.dt
        point['Dt']=self.param['Dt']
        point['eps']=self.param['eps']
        point['dx']=self.dx
        point['mesh']=self.mesh
        point.append()
        table.flush()
        h5file.close()
    



 
