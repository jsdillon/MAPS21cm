# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class GlobalSkyModel:
    """Computes GSM from 3 principal components appropriately weighted. 
    
    Works exactly as the fortran code is described in the comments. I couldn't directly check because I couldn't get the fortran to compile. TODO later.

    Takes the frequency (in MHz), the location of the HEALPIX .fits files (which are in my github under ObservationData/GSM), and the HealPIX NSIDE desired.
    
    Parameters
    ---------------
    freq: float
        Frequency in MHz.
    GSMlocation: str
        Path to GSM .fits and components.dat files.
    GSMNSIDE: int
        Power of 2 between 8 and 512.
        
    Class Members
    ----------
    hpMap: numpy array
        Healpix map of the GSM in galactic coordinates
    """
    
    def __init__(self,freq,GSMlocation,GSMNSIDE):
        self.freq = freq
        self.NSIDE = GSMNSIDE        
        GSMComponents = np.asarray([hp.read_map(GSMlocation + "gsm" + str(i+1) + ".fits" + str(GSMNSIDE),verbose=False) for i in range(3)])
        components = np.loadtxt(GSMlocation + "components.dat")
        temperature = 10**(interp1d(np.log10(components[:,0]), np.log10(components[:,1]), kind='cubic')(np.log10(freq))) #cubic spline interpolation in log(T)
        weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='cubic')(freq) for i in range(3)]) #cubic spline interpolation for weights
        self.hpMap = temperature*np.dot(weights,GSMComponents)