import Solver

class LinearSolver(Solver.Solver):

	def __init__(self,parameters=None):
		Solver.Solver.__init__(self,parameters)

	def solve(self,rhs):
	        """
               	Solve a linear system for the right-hand side value rhs.  
		The method returns (solution,status)
		Status means:
			0: method has converged
			k>0: method has not converged for some reason
                """
		raise NotImplementedError, "a solve() method needs " + \
			"to be implemented for each linear solver"
	
	def new_build(self):
	    pass
