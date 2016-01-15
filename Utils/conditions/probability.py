import scipy

def condition(point,art_index):

    l = len(point.param['free'])
    a = len(point.param['artificial'])
    
    ones = scipy.ones((len(point.u),))
    
    result = {}
    result['column'] = ones
    result['row'] = ones
    result['d'] = scipy.zeros((l+a,),scipy.float64)
    result['eq_term'] = ones*point.lambd[art_index]
    result['res'] = 0

    return result
        
