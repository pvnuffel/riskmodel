class Solver(object):

    def __init__(self,parameters=None):
        if parameters==None:
            self.param = self.getDefaultParameters()
        else:
            self.param=parameters

    def getDefaultParameters():
        raise NotImplementedError, \
            "The default parameters have to be set for each solver"
    getDefaultParameters = staticmethod(getDefaultParameters)
    
    def setPoint(self,point):
        self.point=point
    
    def new_build(self):
        pass

