import numpy as np

def movingaverage1D(data, n_pts):
    ''' calculate a moving average
    data = 1D array to average
    n_pts = number of points to average over'''
    
    #take moving average
    window = np.ones(int(n_pts))/float(n_pts)
    result = np.convolve(data, window, 'same')
    
    #take care of points at ends of trajectories
    result[:n_pts/2] = np.mean(data[:n_pts/2])
    result[-n_pts/2:] = np.mean(data[-n_pts/2:])
        
    return result
