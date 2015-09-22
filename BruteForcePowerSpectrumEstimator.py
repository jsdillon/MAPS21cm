import numpy as np
import healpy as hp
import math
import os
import sys
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
        self.s, self.times, self.ps, self.Dmatrix, self.PSF, self.coaddedMap, self.pointSourcePSF = MapMats.loadAllResults(resultsDirectory)

def loadAllMapmakingResults(s,freqs):
    allMapmakingResults = []    
    for freq in freqs:
        print "Now loading " + s.resultsFolderFormat + "{:.3f}".format(freq) + "/"
        allMapmakingResults.append(MapmakingResult(s.resultsFolderFormat + "{:.3f}".format(freq) + "/"))
    return allMapmakingResults



def ComputeNoiseCovarianceInverse(s,allMapmakingResults):
    print "testing"
    






def BruteForcePowerSpectrumEstimator(configFile = "configuration.txt", **kwargs):
    s = Specifications(os.path.dirname(os.path.abspath(__file__)), "/" + configFile)    
    s.OverrideSpecifications(kwargs)
    coords = Geometry.Coordinates(s)    
    freqs = s.frequencyList[(s.frequencyList >= s.frequencyRange[0]) * (s.frequencyList <= s.frequencyRange[1])]
    coords.computeCubeCoordinates(s,freqs)    
    
    plt.scatter(180/np.pi*coords.pixelRAs[coords.mapIndices], 180/np.pi*coords.pixelDecs[coords.mapIndices])
    plt.figure()    
    plt.scatter(180/np.pi*coords.pixelRAs[coords.mapIndices], coords.xCoords[0,:])
    plt.figure()    
    plt.scatter(180/np.pi*coords.pixelDecs[coords.mapIndices], coords.yCoords[0,:])
    #allMapmakingResults = loadAllMapmakingResults(s,freqs)
    
#    facetCenterIndex = coords.extendedIndices[coords.extendedIndexOfFacetCenter]
 #   for index in np.arange(facetCenterIndex - 5, facetCenterIndex + 5):
 #       print str(index) + ": " + str(180/np.pi*coords.pixelRAs[index]) + ": " + str(180/np.pi*coords.pixelDecs[index])
        
    
    #Compute C
    #Compute Cinv
    #Compute 
    #Compute Q_alpha
    #Compute phat
    #Compte Fisher
    
        
        
    
    
    
    
    

if __name__ == "__main__":
    BruteForcePowerSpectrumEstimator()