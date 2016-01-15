import scipy
import tables

import Point.Point as Point

class Continuer(object):
    
    def __init__(self,points=None,params=None):
        """
        input:
        ======
            points: a list of points to start the branch
                    the minimum is 2
                    note that the points should have sufficient free
                    parameters
            params: dict
                    look at the docstring of getDefaultParameters()
                    to find out which parameters there are

        behaviour:
        ==========
            This class implements a branch structure, containing a list
            of points obtained through continuation.

        """
        for i in range(len(points)):
            points[i].param['free']+=params['free']   
            points[i].param['extra_condition']+=\
                    params['extra_condition']    
            points[i].extra_condition+=params['extra_condition']    
        self.points = points
        self.bound_condition = params['bound_condition']
        if params == None:
            params = self.getDefaultParameters()
        self.param = params

    def getDefaultParameters():
        """
        Returns a dictionary of default parameters
            
            plotx : measure to plot along x axis on bifurcation diagram
                    default = lambd[0]
            ploty : measure to plot along y axis on bifurcation diagram
                    default = norm(u)
                Note: a measure is a dictionary, containing a field 
                and a func, to indicate what should be plotted
            
            growth_factor:
                    indicates the growth for the pseudo-arclength step
                    length after each successful step
                    default = 1.2
            print : *'short'* -- a line per corrected point
                    'none'  -- no printing 
        """    
        param = {}
        
        plotx = {}
        plotx['field'] = 'lambd'
        plotx['func'] = 0
        param['plotx']=plotx
        ploty = {}
        ploty['field'] = 'u'
        ploty['func'] = 'norm'
        param['ploty']=ploty

        param['growth_factor']=1.2
        param['print']='short'
        param['bound_condition']=[]
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)

    def addPoint(self,point):
        """ adds point to the branch """
        self.points.append(point)

    def removeLast(self):
        """ removes the last point from the branch"""
        self.points = self.points[:-1]

    def bcontinue(self,nb_points):
        """ 
        performs a continuation for this branch, adding nb_points
        points (or less, if there are failed Newton iterations)
        """
        # internal parameters        
        tries = 0
        fail = 0
        rjct = 0
        successful = True
        kontinue = True
        bound = False

        growth = self.param['growth_factor']

        l = len(self.points)
        if l <= 1:
            raise ValueError, "There should be two points in the branch"
        while kontinue and tries <= nb_points:
            tries +=1
            l = len(self.points)
            prev_point = self.points[l-2]
            last_point = self.points[l-1]

            secant = last_point.axpy(-1,prev_point)
            dist = scipy.linalg.norm(secant.u)+\
                    scipy.linalg.norm(secant.lambd)
            if successful and not bound:  
                # het vorige punt was successvol berekend
                steplength = growth*dist
            elif not successful:
                steplength = -dist/2
            elif bound:
                # we hebben de rand geraakt
                # nu moeten we de stap verkleinen en zien of het zo gaat
                steplength = steplength/2
                
            new_point = secant.axpy(-steplength/dist,last_point)
            new_point.secant = {}
            new_point.secant['u']=secant.u
            new_point.secant['lambd']=secant.lambd
            status = new_point.correct()
            if status == 0:
                new_success = True
            else: 
                new_success = False
            # check if the bound is reached
            bound = self.boundCrossed(new_point)
            
            if not self.param['print']=='none':
                if bound:
                    bound_hit = "| BOUNDARY HIT!"
                else:
                    bound_hit = ""
                print "Continuation | status : ", status,\
                    " u : " , scipy.linalg.norm(new_point.u),\
                    " lambd: ", new_point.lambd, bound_hit

            if not new_success :
                fail += 1
                if not successful:  # this means that we had two failures
                                    # in a row
                    kontinue = False
                    rjct = rjct + 1
            else :  # if new_success
                if successful and not bound:
                    self.addPoint(new_point)
                elif not successful:   # vorige was geen success 
                        # dus we zijn tussen de twee laatste bezig
                    self.removeLast()
                    self.addPoint(new_point)
                    self.addPoint(last_point)
            successful = new_success
        succ = tries - fail
        return succ,fail,rjct

    def boundCrossed(self,point):
        bound = False
        for bc in self.bound_condition:
            bound = bound or bc(point)
            if bound : 
                break
        return bound
               
               
    def getPlot(self):
        """ 
        Prints a line containing the x and y coordinates
        (corresponding to plotx and ploty in self.param) to plot
        in a bifurcation diagram.
        Also returns these x and y coordinates as a vector.
        """
        x = scipy.zeros((len(self.points),),scipy.float64)
        y = scipy.zeros((len(self.points),),scipy.float64)
        for i in range(len(self.points)):
            x[i] = self.get(self.points[i],"x")
            y[i] = self.get(self.points[i],"y")
            print x[i],y[i]
        return x,y

    def get(self,point,option):
        if option == "x":
            plt = self.param['plotx']
        elif option == "y":
            plt = self.param['ploty']
        else: 
            raise ValueError, "You gave a "\
                + "nonsense option to get()!!"

        field=point.__getattribute__(plt['field'])
        func=plt['func']
        if type(func)==int:
            return field[func]
        elif func == "norm":
            return scipy.linalg.norm(field)

    def save(self,file,group,table):
        for p in self.points:
            p.save(file,group,table)

def branch_read(filename,groupname,tablename):
    h5file = tables.openFile(filename,"r")
    exec "table = h5file.root." + groupname + "." + tablename
    points = []
    for row in table:
        point = {}
        for name in table.colnames:
            point[name]=row[name]
        points.append(point)
    h5file.close()
    return points
            
        


