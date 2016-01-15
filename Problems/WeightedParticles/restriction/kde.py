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
        edges = scipy.r_[domain[0],(grid[1:]+grid[:-1])/2.,domain[-1]]
        # estimating the cumulative density to make the estimation conservative
        rho = statistics.pdf(x,edges[1:-1],weight=weight,h=self.h)
        return rho
    
