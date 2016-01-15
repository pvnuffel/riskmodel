import scipy

class Histogram(object):
    
    def __init__(self,param=None):
        pass
    
    def getDefaultParameters():
        param = {}
        return param
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def restrict(self,x,grid,domain,w=None):
        edges = scipy.r_[domain[0],(grid[1:]+grid[:-1])/2.,domain[-1]]
        rho,edges = scipy.histogram(x,edges,weights=w,normed=True)
        return rho
    
    def getBins(self,x,grid,domain):
        edges = scipy.r_[domain[0],(grid[1:]+grid[:-1])/2.,domain[-1]]
        bins = scipy.digitize(x,edges)-1
        return bins
