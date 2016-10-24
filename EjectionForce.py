''' This is an implementation of the calculation outlined in Purohit, Kindev, Phillips, 
    Mechanics of DNA packaging in viruses, PNAS, vol 100, no 6, 3173-3178, 2003. http://www.pnas.org/content/100/6/3173
    The authors state (just below equation 19) that the force using a cylindrical geometry is essentially 
    the same as if a sphere with the same volume is used. See Figure 2 for an comparison of ds for each geometry.
    
    L is the contour length of the genome in the capsid.
    F0 = characteristic force (actually a pressure) of DNA strand-strand repulsion. 
    c = characteristic length scale of DNA strand-strand repulsion. 
    Rout = radius of capsid cylinder or sphere.
    z = height of cylindrical capsid (not used for sphere calculation).
    psi_p = DNA persistance length. 
    kBT = Thermal energy. 
    geom = capsid geometry. Cylindrical 'cyl', or spherical 'sph'.
    
    You should be able to use any self-consistent unit set for the parameters.
'''

import numpy as np
from scipy.optimize import brentq

def ejection_force_calc(L, F0, c, Rout = 19.4, z = 37.9, psi_p = 50., kBT = 1.38e-23*(273.15+37) *1e12*1e9, geom = 'cyl'):
    ''' Ejection force calculation.
    All input parameters except for geom should be floats. geom should be a string.
    Returns ejection force. '''
    
    if L <= 0: # L <= 0 causes the zero-finder below to fail, but I can do this force calculation by hand.
        return 0.0
    
    if geom == 'cyl':
        def R(ds): # this is the expression for R that is in between equation 6 and 7. 
            return Rout*np.sqrt(1 - np.sqrt(3)*ds**2*L/(2*np.pi*z*Rout**2) )
                
        def to_zero(ds): #this is left_hand_side - right_hand_size of equation 17
            return np.sqrt(3)*F0*np.exp(-ds/c) - ( psi_p*kBT/R(ds)**2/ds**2 - 2*psi_p*kBT/ds**2*np.log(Rout/R(ds))/(Rout**2-R(ds)**2) )
            
        # Prevent the minimizer from taking sqrt of a negatve number when calculating R(ds).
        # The below expression was determined by finding when the argument of the sqrt in R(ds) becomes negative. 
        # This equivlaent to keeping R(ds) > 0.
        ds_max = 0.999999999*Rout*np.sqrt(2*np.pi*z/np.sqrt(3)/L)
        
    elif geom == 'sph':
        def R(ds): # Determined by solving L(R) for a sphere (given in Table 1) for R.
            return Rout*np.sqrt(1 - (3*np.sqrt(3)*L*ds**2/8/np.pi)**(2./3.)/Rout**2 )
                
        def to_zero(ds): #this is left_hand_side - right_hand_size of equation 18
            return np.sqrt(3)*F0*np.exp(-ds/c) - ( psi_p*kBT/R(ds)**2/ds**2 + 3*psi_p*kBT/ds**2 * ( 1/(Rout**2-R(ds)**2) + Rout/(Rout**2-R(ds)**2)**1.5*np.log((Rout-np.sqrt(Rout**2-R(ds)**2))/R(ds)) ) )
            
        # Prevent the minimizer from taking sqrt of a negatve number when calculating R(ds).
        ds_max = 0.999999999*np.sqrt(8*np.pi*Rout**3/3/np.sqrt(3)/L)     
    
    # Note that above there are instances where Rout**2-R(ds)**2 appears in the denominator.
    # However, as long as ds > 0, then this expression will be  > 0 so we do not have to worry about dividing
    # by zero as long as ds_min > 0.
    ds_min = Rout*1e-4 #make ds_min much smaller than capsid radius (I find anything in the range 0.0002 nm - 0.1 nm works)
    
    # calculate ds by finding a zero
    # returns the zero location in the interval [ds_min,ds_max].
    ds = brentq(to_zero, ds_min, ds_max)
    
    F = np.sqrt(3)*F0*(c**2 + c*ds)*np.exp(-ds/c)+psi_p*kBT/2/R(ds)**2 #equation 19
    
    return F


def ejection_force(L, F0, c, Rout = 19.4, z = 37.9, psi_p = 50., kBT = 1.38e-23*(273.15+37) *1e12*1e9, geom = 'cyl'):
    '''wrapper function for calculation to deal with arrays of L. 
    TODO: generalize this wrapper so any number of input parameters could be arrays. '''
    
    
    if isinstance(L,np.ndarray): # if L is an array
        #explicitly loop through each element of the array since the brentq minimizer doesn't work on arrays.
        result = np.empty(len(L))
        for i in range(len(L)):
            result[i] = ejection_force_calc(L[i], F0, c, Rout = Rout, z = z, psi_p = psi_p, kBT = kBT, geom = geom)
    else: #if L is a scalar
        result = ejection_force_calc(L, F0, c, Rout = Rout, z = z, psi_p = psi_p, kBT = kBT, geom = geom)
    return result
 

def test_force_calculations():
    '''test calculation of force to make sure it works.
    This should reproduce the theory curves in Figure 3 of the paper for phage phi29.
    It also plots the force using a spherical geometry, assuming the volume of the capsid is conserved.
    The force for each geometry should be very similar, with the force in the sphere being slightly higher.
    '''
    
    import matplotlib.pyplot as plt
    
    #in units of nm an pN
    F0s = np.array([2,4.1,6])*55000. #pN/nm^2
    Rout = 19.4 # nm
    z = 37.9 # nm
    L_max = 6580. # nm
    psi_p = 50. # nm
    c = 0.27 # nm
    kBT = 1.38e-23*(273.15+37) *1e12*1e9 # convert to pN and nm
    
    Rout_sph = (3.*Rout**2*z/4.)**(1./3) #radius of sphere that has the same volume as the cylinder

    L = np.linspace(0,L_max, 100)

    plt.figure()
    for F0 in F0s:
        F_cyl  = ejection_force(L, F0, c, Rout = Rout, z = z, psi_p = psi_p, kBT = kBT, geom = 'cyl')
        F_sph = ejection_force(L, F0, c, Rout = Rout_sph, psi_p = psi_p, kBT = kBT, geom = 'sph')

        plt.plot(L/L_max*100, F_cyl, label = 'cyl. F0: ' + str(F0/55000))
        plt.plot(L/L_max*100, F_sph, label = 'sph. F0: ' + str(F0/55000))
    plt.ylim([-10,90])
    plt.xlim([0,120])
    plt.ylabel('Force (pN)')
    plt.xlabel('Percent packed')
    plt.legend(loc = 'upper left')
    
