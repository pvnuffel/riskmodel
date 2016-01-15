import scipy
import tables

class System(object):

    def __init__(self,u,lambd,parameters=None):
        """
        initializes a system with the current state and parameter values
        """
        if parameters == None:
            self.param = self.getDefaultParameters()
        else:
            self.param = parameters
        self.setState(u,lambd)
    
    def getDefaultParameters():
        return {}
    getDefaultParameters = staticmethod(getDefaultParameters)
    
    def setState(self,u,lambd):
        """
        This method is used to set the current system state
        in a consistent manner
        """
        self.u = u
        self.lambd = lambd
    
    def getResidual(self):
        """
        This method determines what the residual for the system is.
        """
        raise NotImplementedError, \
            "System.getResidual() needs to be implemented in " +\
            "subclass"
    
    def computeJacobian(self):
        """
        returns the Jacobian of the system for the given state and
        parameter vectors
        """
        raise NotImplementedError, "System.computeJacobian() should be "+\
            "implemented in the concrete (equation-based) system"
    
    def applyJacobian(self,v):
        """
        Return the Jacobian-vector product of the
        system Jacobian with the vector v
        """
        raise NotImplementedError, "System.applyJacobian() should be" +\
            "implemented for the concrete type of systems"
    
    def solvePreconditioner(self,rhs,rows,cols):
        """ 
        solves a linear system with the preconditioning matrix
        
        input:
        =====
            rhs 	contains the right-hand side of the system to solve
            rows	contains a number of extra rows that come from external
                    constraints
            cols 	contains a number of extra columns that contain entries 
                    stemming from free parameters
        """
        raise NotImplementedError, "System.solvePreconditioner() should "+\
            "be implemented for a concrete system"
    
    def applyPreconditioner(self,v):
        """ 
        solves a linear system with the preconditioning matrix
        
        input:
        =====
            rhs 	contains the right-hand side of the system to solve
            rows	contains a number of extra rows that come from external
                    constraints
            cols 	contains a number of extra columns that contain entries 
                    stemming from free parameters
        """
        raise NotImplementedError, "should be implemented for a concrete system"
    
    def getParameterDerivative(self,i):
        """
        This method returns the derivative with respect to parameter i
        """
        raise NotImplementedError, \
            "System.getParameterDerivative() needs to be " +\
            "implemented by a subclass"
    
    def getHDF5Description(self):
        u_shape = scipy.shape(self.u)
        lambd_shape = scipy.shape(self.lambd)
        class SystemSave(tables.IsDescription):
            u = tables.FloatCol(shape = u_shape)
            lambd = tables.FloatCol(shape = lambd_shape)
        return SystemSave
        
    
    def save(self,file,group,table):
        h5file = tables.openFile(file,mode="a")
        table = self.getTable(h5file,group,table)
        point=table.row
        point['u']=self.u
        point['lambd']=self.lambd
        print point['lambd']
        point.append()
        table.flush()
        h5file.close()
    
    def getTable(self,h5file,group,table):
        try:
            exec "table = h5file.root." + group + "." + table
        except tables.exceptions.NoSuchNodeError:
            group = h5file.createGroup("/",group,group)
            system = self.getHDF5Description()
            table = h5file.createTable(group,table,system,"system")
        return table
    
                
