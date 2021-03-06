import scipy

class Sampler(object):    
    def __init__(self,param=None):
        if param == None :
            param = Sampler.getDefaultParameters()
        self.param = param
        self.rand = scipy.random.mtrand.RandomState()
        self.seed(param['seed'])
    
    def getDefaultParameters():
        param = {}
        param['seed']=0
        return param
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def lift(self,rho,grid,N):
        cum,edges = self.rho2cum(rho,grid)
        y = self.rand.uniform(size=N)
        print "sampling density with seed : ", self.cur_seed, "first sample: ", y[0]
        return self.inv_cum(y,cum,edges)
    
    def rho2cum(self,rho,grid):
        ngrid = len(grid)
        cum = scipy.zeros((ngrid+1,))
        edges = scipy.zeros_like(cum)
        # grid values are box centers
        edges[0]=3./2.*grid[0]-1./2.*grid[1]
        edges[1:ngrid]=(grid[:ngrid-1]+grid[1:ngrid])/2.
        edges[-1]=3./2.*grid[-1]-1./2.*grid[-2]
        
        cum[0]=0.
        for n in range(ngrid):
            if rho[n]<0 :
                print "encountered a negative density : ", rho[n]
            cum[n+1]=cum[n]+rho[n]*(edges[n+1]-edges[n])
        print "estimating cumulative density : last value ", cum[-1], " (should be 1.)"
        return cum,edges
    
    def inv_cum(self,y,cum,edges):
        x = scipy.zeros_like(y)
        for i in range(len(y)):
            x[i]=self.inv_cum_single_particle(y[i],cum,edges)
        return x
    
    def inv_cum_single_particle(self,yn,cum,edges):
        i = 0
        while (yn>cum[i]):
            i+=1
        m = (cum[i]-cum[i-1])/(edges[i]-edges[i-1])
        xn = (yn-cum[i-1])/m+edges[i-1]
        return xn
    
    def seed(self,s):
        print "seeding the density sampler with seed : ", s
        self.cur_seed = s
        self.rand.seed(self.cur_seed)
    

# Here we put a number of tests :
if __name__ == "__main__":
    
    # construction of initial grid and density
    xL = 0.
    xR = 20.
    Dx = 0.2
    grid = scipy.arange(xL+Dx/2.,xR,Dx)
    ngrid = len(grid)
    rho = scipy.ones_like(grid)/(xR-xL)

    # construction of sampler 
    param = Sampler.getDefaultParameters()
    param['N']=20
    sampler = Sampler(param)

    # construction of cumulative distribution (should be 1)
    cum,edges = sampler.rho2cum(rho,grid)
    print "test cumulative: ", cum[-1],edges[0],edges[-1]

    print cum
    print edges
    # construction of positions, since cumulative is uniform, we should 
    # get back the random numbers that we put in
    y = scipy.random.uniform(size=sampler.param['N'])
    x = sampler.inv_cum(y,cum,edges)
    
    # needs to be a rescaled version of y, rescaled by the size of the domain
    print "ratio: ", x/y
    
    # check if the result changes with changing order of x
    ysort = scipy.sort(y)
    indices = scipy.argsort(y)
    xsort = sampler.inv_cum(ysort,cum,edges)
    print scipy.amax(xsort-x[indices])   

    
