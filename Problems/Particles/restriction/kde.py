
# import sys 
# sys.path.append('../')
# sys.path.append('../../')
# sys.path.append('../../../')
# sys.path.append('../Modules/')
import statistics
import scipy

class KDE(object):
    
    def __init__(self,param):
        self.h = param['h']
    
    def getDefaultParameters():
        param = {}
        param['h']=1e-1
        return param
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def restrict(self,x,grid,domain):
        edges = scipy.r_[domain[0],(grid[1:]+grid[:-1])/2.,domain[-1]]   #concatenate
        # estimating the cumulative density to make the estimation conservative
        cum = statistics.cpdf(x,edges[1:-1],h=self.h)
        rho = scipy.zeros_like(grid)
        rho[0]=cum[0]/(edges[1]-edges[0])
        rho[1:-1]=(cum[1:]-cum[:-1])/(edges[2:-1]-edges[1:-2])
        rho[-1]=(1.-cum[-1])/(edges[-1]-edges[-2]) 
        return rho
    
