import scipy

class Histogram(object):
    
    def __init__(self):
        pass
    
    def restrict(self,x,grid,domain):
        edges = scipy.r_[domain[0],(grid[1:]+grid[:-1])/2.,domain[-1]]
        rho,edges = scipy.histogram(x,edges,density=True)
        return rho
    

