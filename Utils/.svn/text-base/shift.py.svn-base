import scipy

def shift_back_front(x,phi,neq,dx):
	n = len(x)/neq
	u = scipy.reshape(x,(neq,n))
	u_left = scipy.transpose([u[:,0]])
	u_right = scipy.transpose([u[:,-1]])
	u_ext = scipy.c_[u_left,u,u_right]
	dudx = 1./(2*dx)*(u_ext[:,2:]-u_ext[:,:-2])
	#u_ext = scipy.c_[u_left,u_left,u,u_right,u_right]
	#dudx = 1./(12*dx)*(-u_ext[:,4:]+8*u_ext[:,3:-1]-8*u_ext[:,1:-3]+u_ext[:,:-4])
	x_shifted=scipy.reshape(u+phi*dudx,(neq*n,))
	return x_shifted

def shift_back_pulse(x,phi,neq,mesh,dx):
	
	L=mesh[-1]-mesh[0]+dx
	u=scipy.reshape(x,(neq,scipy.shape(mesh)[0]))
	uhat=1./len(mesh)*scipy.fftpack.rfft(u)
	u_shifted=scipy.outerproduct(uhat[:,0],\
		scipy.ones(scipy.shape(mesh),\
				scipy.float64))
	if len(mesh)%2==0:
		u_shifted=u_shifted+scipy.outerproduct(uhat[:,-1],\
			(-1.)**scipy.arange(len(mesh)))
	for k in range(1,(len(mesh)-1)/2+1):
	  u_shifted=u_shifted+\
	    2*scipy.outerproduct(\
		    uhat[:,2*k-1],scipy.cos(\
		2*k*scipy.pi/L*(mesh-mesh[0]+phi)))-\
	    2*scipy.outerproduct(\
		    uhat[:,2*k],scipy.sin(\
		2*k*scipy.pi/L*(mesh-mesh[0]+phi)))
	x_shifted=scipy.reshape(u_shifted,\
		(neq*scipy.shape(mesh)[0],))
	return x_shifted


