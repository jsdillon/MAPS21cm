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
from scipy.integrate import dblquad
from scipy.special import jv
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



def computeNoiseCovariance(allMapmakingResults,coords):
    """This function computes the noise covariance. The indices are CN[voxel_i, voxel_j]."""
    CN = np.zeros((coords.nVoxels,coords.nVoxels))    
    for channel in range(coords.nFreqs):
        CNofFreq = np.dot(allMapmakingResults[channel].PSF[:,coords.mapIndexLocationsInExtendedIndexList],np.transpose(np.diag(allMapmakingResults[channel].Dmatrix)))
        CN[channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels, channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels] = CNofFreq
    return CN
    

def computePtransCinvP(allMapmakingResults,Cinv,coords):
    """This funciton computes Ptranspose * Cinv * P where P is approximated to be square (essentially PSFextension factor = 1)"""
    PtransCinvP = np.zeros((coords.nVoxels,coords.nVoxels))    
    for channel in range(coords.nFreqs):
        PSF = allMapmakingResults[channel].PSF[:,coords.mapIndexLocationsInExtendedIndexList]
        CinvHere = Cinv[channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels, channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels]
        PtransCinvPofFreq = np.transpose(PSF).dot(CinvHere).dot(PSF)
        PtransCinvP[channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels, channel*coords.nFacetPixels:(channel+1)*coords.nFacetPixels] = PtransCinvPofFreq
    return PtransCinvP

def computeFourierBinEdges(coords):
    """This function retruns the edges of kParaBins/kPerpBins with indices [binNumber][start/stop]."""
    perpLength = ((np.max(coords.xCoords) - np.min(coords.xCoords))**2 + (np.max(coords.yCoords) - np.min(coords.yCoords))**2)**.5
    paraLength = np.max(coords.losCoords) - np.min(coords.losCoords)
    deltaKPara = 2*np.pi / paraLength
    deltaKPerp = 2*np.pi / perpLength
    kParaBinEdges = [(deltaKPara*n, deltaKPara*(n+1)) for n in range(int(np.ceil(len(coords.freqs)/2.0)))]
    kPerpBinEdges = [(deltaKPerp*n, deltaKPerp*(n+1)) for n in range(int(np.ceil(coords.nFacetPixels**.5 / 2)))]
    return kParaBinEdges, kPerpBinEdges 

def calculateQMatrices(coords, kParaBinEdges, kPerpBinEdges):
    """Computes Q_alpha matrix assuming delta function pixelization. Indices are Q[kParaBin, kPerpBin, voxel_i, voxel_j].""" 
    QMatrices = np.zeros((len(kParaBinEdges), len(kPerpBinEdges), coords.nVoxels, coords.nVoxels))
    rPara = np.zeros((coords.nVoxels, coords.nVoxels))
    rPerp = np.zeros((coords.nVoxels, coords.nVoxels))
    for iLos in range(coords.nFreqs):
        for iPix in range(coords.nFacetPixels):
            i = iLos*coords.nFacetPixels + iPix
            for jLos in range(coords.nFreqs):
                for jPix in range(coords.nFacetPixels):
                    j = jLos*coords.nFacetPixels + jPix
                    rPara[i,j] = (coords.losCoords[iLos,iPix] - coords.losCoords[jLos,jPix])
                    rPerp[i,j] = ((coords.xCoords[iLos,iPix] - coords.xCoords[jLos,jPix])**2 + (coords.yCoords[iLos,iPix] - coords.yCoords[jLos,jPix])**2)**.5

    for kParaBin in range(len(kParaBinEdges)):    
        for kPerpBin in range(len(kPerpBinEdges)):
            paraTerm = (2*np.pi**2)**-1 * rPara**-1 * (np.sin(kParaBinEdges[kParaBin][1] * rPara) - np.sin(kParaBinEdges[kParaBin][0] * rPara))
            paraTerm[rPara == 0] = (2*np.pi**2)**-1 * (kParaBinEdges[kParaBin][1] - kParaBinEdges[kParaBin][0])            
            perpTerm = rPerp**-2 * (kPerpBinEdges[kPerpBin][1] * rPerp * jv(1, kPerpBinEdges[kPerpBin][1] * rPerp) - kPerpBinEdges[kPerpBin][0] * rPerp * jv(1, kPerpBinEdges[kPerpBin][0] * rPerp))
            perpTerm[rPerp == 0] = (kPerpBinEdges[kPerpBin][1]**2 - kPerpBinEdges[kPerpBin][0]**2) / 2.0
            QMatrices[kParaBin,kPerpBin] = paraTerm * perpTerm
            #TODO: consider series expansion for very small r

    return QMatrices 


def calculateFisherMatrix(PtransCinvP, QMatrices):
    """Calculates the Fisher Matrix as .5*tr[Q_alpha * Cinv * Q_beta * Cinv]. Indices loop over kPerp bins fastest and kPara bins slower.""" 
    kParaBins = QMatrices.shape[0]
    kPerpBins = QMatrices.shape[1]
    Fisher = np.zeros((kParaBins*kPerpBins, kParaBins*kPerpBins))
    for kParaBinAlpha in range(kParaBins):
        for kPerpBinAlpha in range(kPerpBins):
            for kParaBinBeta in range(kParaBins):
                for kPerpBinBeta in range(kPerpBins):
                    traceAnswer = .5 * np.trace(QMatrices[kParaBinAlpha][kPerpBinAlpha].dot(PtransCinvP).dot(QMatrices[kParaBinBeta][kPerpBinBeta]).dot(PtransCinvP))
                    Fisher[kParaBinAlpha*kPerpBins + kPerpBinAlpha][kParaBinBeta*kPerpBins + kPerpBinBeta] = traceAnswer            
    return Fisher

def BruteForcePowerSpectrumEstimator(configFile = "configuration.txt", **kwargs):
    print "stuff"
    #Q_alpha
    
    

    #Compute C
    #Compute Cinv
    #Compute 
    #Compute Q_alpha
    #Compute phat
    #Compte Fisher
    
        
        
    
    
        #plt.plot(freqs,coords.comovingDistances,'.')
    
#    plt.scatter(180/np.pi*coords.pixelRAs[coords.mapIndices], 180/np.pi*coords.pixelDecs[coords.mapIndices])
#    plt.figure()    
#    plt.scatter(180/np.pi*coords.pixelRAs[coords.mapIndices], coords.xCoords[0,:])
#    plt.figure()    
#    plt.scatter(180/np.pi*coords.pixelDecs[coords.mapIndices], coords.yCoords[0,:])
#
#

#    
    
#    facetCenterIndex = coords.extendedIndices[coords.extendedIndexOfFacetCenter]
#    for index in np.arange(facetCenterIndex - 5, facetCenterIndex + 5):
#        print str(index) + ": " + str(180/np.pi*coords.pixelRAs[index]) + ": " + str(180/np.pi*coords.pixelDecs[index])
            
    
    
    

if __name__ == "__main__":
    BruteForcePowerSpectrumEstimator()

configFile = "configuration.txt"
kwargs = {}

s = Specifications(os.path.dirname(os.path.abspath(__file__)), "/" + configFile)    
s.OverrideSpecifications(kwargs)
coords = Geometry.Coordinates(s)    
freqs = s.frequencyList[(s.frequencyList >= s.frequencyRange[0]) * (s.frequencyList <= s.frequencyRange[1])]
coords.computeCubeCoordinates(s,freqs)

allMapmakingResults = loadAllMapmakingResults(s,freqs)
CN = computeNoiseCovariance(allMapmakingResults,coords)
Cinv = np.linalg.inv(CN)
PtransCinvP = computePtransCinvP(allMapmakingResults,Cinv,coords)
#PtransCinvP = np.identity(len(PtransCinvP))

kParaBinEdges, kPerpBinEdges = computeFourierBinEdges(coords)
QMatrices = calculateQMatrices(coords, kParaBinEdges, kPerpBinEdges)
Fisher = calculateFisherMatrix(PtransCinvP, QMatrices)





#%% 
def pColorKParaKPerp(whatToPlot,kParaBinEdges,kPerpBinEdges,title):
    fig = plt.figure()
    kParaForPColor = np.append(np.asarray(kParaBinEdges)[:,0],(kParaBinEdges[-1][1]))
    kPerpForPColor = np.append(np.asarray(kPerpBinEdges)[:,0],(kPerpBinEdges[-1][1]))
    plt.pcolor(kPerpForPColor, kParaForPColor, whatToPlot, norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.xlabel('$k_\perp$ Mpc$^{-1}$')
    plt.ylabel('$k_\|$ Mpc$^{-1}$')
    plt.title(title)
    return fig


#verticalErrors = np.reshape(np.transpose(np.diag(np.linalg.inv(Fisher))),(len(kPerpBinEdges),len(kParaBinEdges)))
FisherDiag = np.reshape(np.diag(Fisher),(len(kParaBinEdges),len(kPerpBinEdges)))
fig = pColorKParaKPerp(FisherDiag,kParaBinEdges,kPerpBinEdges,'Fisher Diagonal')

