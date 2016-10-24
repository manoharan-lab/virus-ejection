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

def RemovePauses(datas, pauses):
    ''' remove pauses from data
    datas = list of 1D arrays to remove pauses from
    pauses = list of lists of the start and stop indecies of each pause. 
        Pauses can be listed in any order'''
	
    dataps = []
    for i in range( 0,len(datas) ):
        data = datas[i]
        datal = list(data)
        pause = pauses[i].sort().reverse()
        for j in range ( 0,len(pause),2 ):
            del datal[pause[j]:pause[j+1]]
        dataps.append(np.array(datal))
	
    return dataps
