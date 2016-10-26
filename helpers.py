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
        pause = map(int,pauses[i])
        pause.sort()
        pause.reverse()
        for j in range ( 0,len(pause),2 ):
            del datal[pause[j+1]:pause[j]]
        dataps.append(np.array(datal))
    return dataps
    
def load_dict(filename, delimiter = ','):
    '''loads a dictionary of parameters from a file. white space is ignored
    parameters could be a number, string, list of numbers, or list of strings.
    numbers with basic + - math charcters are accepted
    '''
    with open(filename) as f:
        ej_content = f.readlines()    
    params = {}
    for item in ej_content:        
        # remove white space
        item = item.replace(" ", "") 
        item = item.replace("\n", "")
        item = item.replace("\t", "")
        if item != '':
            
            #split into name and value
            name, value = item.split('=')
            
            #split into a list if needed
            if delimiter in value:
                value = value.split(delimiter)
                try: #convert to a list of numbers if possible
                    value = map(float, value)
                except:
                    #see if this is a list of strings or a list of numbers with +-
                    found_az = False
                    for item in value:  #if there aren't a-z or A-Z characters, assume it is numbers with +-
                        if re.search('[a-zA-Z]', item) != None: found_az=True
                    if not found_az: value = str_to_float_plusminus(value)
                                   
            else: #if not a list convert to a number if possible
                if value == 'None': value = None
                try: 
                    value = float(value)
                except: ''
                                      
            params[name] = value
    return params                
        
def str_to_float_plusminus(data):
    ''' converts a string or list of strings with + and - signs to a float or list of floats
    '''
    
    #treat a string input as a list of length 1
    if not isinstance(data,list):
        data = [data]
    
    #go through each item
    for i in range(len(data)):
        
        #parse based on plus sign
        numsplus = data[i].split('+')        
        for j in range(len(numsplus)):
        
            #parse based on minus sign
            numsminus = numsplus[j].split('-')
            
            if len(numsminus) == 1: #if there were no minus signs
                numsplus[j] = float(numsminus[0])
            else:
                if numsminus[0] == '':# if first character is a minus sign 
                    numsplus[j] = -1*float(numsminus[1])
                    k = 2
                else:# if first character isn't a minus sign
                    numsplus[j] = float(numsminus[0])
                    k = 1
                while k < len(numsminus): #subtract all numbers
                    numsplus[j] = numsplus[j] - float(numsminus[k])
                    k = k+1
         
        data[i] = sum(numsplus)  #add up numbers
               
    if len(data) == 1:
        data = data[0]
        
    return data
    
def save_dict(filename, params, delimiter = ','):
    ''' saves a dictionary in text file
    
    params is the dictionary. Each value of params can be a number, string, or list of numbers or strings.
    Each item of each value is converted to a string before being written
    '''

    #open file
    f = open(filename,'w')   
    #go through each key
    for key, value in params.iteritems():
        outline = key + ' = '
        #add each item if value is a list
        if isinstance(value, list) or isinstance(value, np.ndarray):
            for i in range(len(value)):
                outline = outline + str(value[i])
                if i != len(value)-1:
                    outline = outline + delimiter+ ' '
        else: #just add the single item
            outline = outline + str(value)
        
        #write line
        f.write(outline+'\n')
    f.close()
    
def to_odd(value):
    '''coverts value to the nearest odd integer, where value is some scalar
    if value is an even integer the ouput is value+1
    '''
    valuemod2 = value%2
    if valuemod2 < 1:
        value = int(value) + 1
    else:
        value = int(value)
    return value
