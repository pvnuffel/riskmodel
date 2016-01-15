import scipy

def condition(point,art_index):

    l = len(point.param['free'])
    a = len(point.param['artificial'])
    
    n = len(point.u)
    result = {}
    result['column'] = scipy.ones((n,))/n
    result['row'] = 2 * scipy.ones((n,))/n
    result['d'] = scipy.zeros((l+a,),scipy.float64)
    result['eq_term'] = scipy.ones((n,))/n*point.lambd[art_index]
    result['res'] = 0

    return result
        
