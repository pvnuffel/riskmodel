import scipy

def condition(point):
    l = len(point.param['free'])
    a = len(point.param['artificial'])

    result = {}
    result['row'] = scipy.zeros(scipy.shape(point.secant['u']),scipy.float64)
    result['d']=scipy.zeros((l+a,),scipy.float64)
    result['d'][l-1] = 1.
    result['res'] = 0
    return result
	    
