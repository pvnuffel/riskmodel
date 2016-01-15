import scipy

def condition(point):

   dx = point.system.mesh[1]-point.system.mesh[0]
   npoints = len(point.system.mesh)
   u = point.system.u[:npoints]
   
   A = scipy.zeros((2,2),scipy.float64)
   b = scipy.zeros((2,),scipy.float64)
   
   # note: the 0.7 is heuristic.   We need to take the tail, but 
   # the numbers shouldn't be too small to make a fit.
   # ("heuristic" then means tried a few values and this one works.)
   N = int(0.7*npoints)
   
   for i in range(2):
       j = N+i
       A[i,0]=u[j]
       A[i,1]= (-u[j+2]+8*u[j+1]-8*u[j-1]+u[j-2])/(12*dx)
     # A[i,1]= (u[j+1]-u[j-1])/(2*dx)
       b[i] = (u[j+1]-2*u[j]+u[j-1])/dx**2
   x = scipy.linalg.solve(A,b)
   # return True betekent dat we op de rand gebotst zijn.
   # dit geeft False terug boven de minimale snelheid
   return 4*x[0]+x[1]**2<0


