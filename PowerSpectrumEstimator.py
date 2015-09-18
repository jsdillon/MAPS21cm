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

class MapmakingResult:
    def __init__(self,resultsDirectory):
        self.s, self.coords, self.times, self.ps, self.Dmatrix, self.PSF, self.coaddedMap, self.mapNoiseCovariance, self.pointSourcePSF = MapMats.loadAllResults(resultsDirectory)

def PowerSpectrumEstimator(configFile = "configuration.txt", **kwargs):
    s = Specifications(os.path.dirname(os.path.abspath(__file__)), "/" + configFile)    
    s.OverrideSpecifications(kwargs)
    cubeFreqs = s.frequencyList[(s.frequencyList >= s.frequencyRange[0]) * (s.frequencyList <= s.frequencyRange[1])]
    #allMapmakingResults = [MapmakingResult(s.resultsFolderFormat + "{:.3f}".format(freq) + "/") for freq in cubeFreqs]
    allMapmakingResults = []    
    for freq in cubeFreqs:
        print "working on " + s.resultsFolderFormat + "{:.3f}".format(freq) + "/"
        allMapmakingResults.append(MapmakingResult(s.resultsFolderFormat + "{:.3f}".format(freq) + "/"))
        
        
    
    
    
    
    
    

if __name__ == "__main__":
    PowerSpectrumEstimator()