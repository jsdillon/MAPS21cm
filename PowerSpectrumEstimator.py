import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
import scipy.constants as const
plt.close('all')

def PowerSpectrumEstimator(configFile = "configuration.txt", **kwargs):
    
    #Load in everything we need, figure out which LSTs to work with
    scriptDirectory = os.path.dirname(os.path.abspath(__file__))
    s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    s.OverrideMapmakingVariables(kwargs)

    
    
    
    
    
    
    

if __name__ == "__main__":
    PowerSpectrumEstimator()