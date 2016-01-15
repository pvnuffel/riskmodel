import scipy

def condition(point,art_index):

    l = len(point.param['free'])
    a = len(point.param['artificial'])
    
    neq = point.neq
    dx = point.system.dx
    nx = len(point.u)/neq

    y0 = scipy.reshape(point.u,(neq,nx))
    
    left = scipy.zeros((neq,1),scipy.float64)
    right = scipy.zeros((neq,1),scipy.float64)
    left[:,0]=y0[:,0]
    right[:,0]=y0[:,-1]
    u=scipy.c_[left,y0,right]
    
    deriv = 1./(2*dx)*scipy.reshape(scipy.transpose(\
        u[:,2:]-u[:,:-2]),(nx*neq,))
    
    result = {}
    result['column'] = deriv
    result['row'] = deriv*dx
    result['d'] = scipy.zeros((l+a,),scipy.float64)
    result['eq_term'] = deriv*point.lambd[art_index]
    result['res'] = 0

    return result
        
