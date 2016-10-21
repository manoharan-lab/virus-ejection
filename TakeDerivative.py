import numpy as np
import helpers

def diffDerivative(data, n_mov_avg):
    ''' takes derivative of data
    data: 1D array to take derivative of
    n_mov_avg width of moving average to use'''
    
    difference = helpers.movingaverage1D( np.diff(helpers.movingaverage1D(data,n_mov_avg)), n_mov_avg)
    
    return difference
