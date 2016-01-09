# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class GlobalSkyModel:
    """Computes GSM from 3 principal components appropriately weighted. 

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
        GSMComponentsDegraded = np.asarray([np.load(GSMlocation + "component_maps_408locked_NSIDE-" + str(GSMNSIDE) + "_Comp-" + str(comp) + ".npy") for comp in range(3)])
        components = np.loadtxt(GSMlocation + "components.dat")
        temperature = np.exp(interp1d(np.log(components[:,0]), np.log(components[:,1]), kind='cubic')(np.log(freq))) #cubic spline interpolation in log(f), log(T)
        weights = np.asarray([interp1d(np.log(components[:,0]), components[:,i+2], kind='cubic')(np.log(freq)) for i in range(3)]) #cubic spline interpolation for log(f), weights
        self.hpMap = temperature*np.dot(weights,GSMComponentsDegraded)