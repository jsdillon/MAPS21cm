import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
import scipy.constants as const

from Source.GlobalSkyModel import GlobalSkyModel
plt.close('all')


def Mapmaker(freq = 150, configFile = "configuration.txt", overridePSFExtension = -1):
    """This function makes maps from visibilities and also calculates the associated map statistics. 
    
    Saves the results to binary (as pickles or numpy arryas) and returns the folder where they are located"""
    
    #Load in everything we need, figure out which LSTs to work with
    scriptDirectory = os.path.dirname(os.path.abspath(__file__))
    s = Specifications(scriptDirectory, "/" + configFile,freq)
    if overridePSFExtension > 0:
        s.PSFextensionBeyondFacetFactor = overridePSFExtension    
    
    times = Geometry.Times(s)
    times.CutOutUnusedLSTsAndGroupIntoSnapshots(s)
    coords = Geometry.Coordinates(s)
    PBs = PrimaryBeams(s)
    ps = PointSourceCatalog(s,PBs,times)
    
    
    #Simulate or load visibilities
    if s.simulateVisibilitiesWithGSM or s.simulateVisibilitiesWithPointSources:
        visibilities = VisibilitySimulator(s,PBs,ps,times,coords)
    else:
        print "Visibility loading functions are not done."    
        #visibilties = LoadVisibilities(s)
    visibilities *= s.convertJyToKFactor
    
    #Prepare visibilities
    Geometry.rephaseVisibilitiesToSnapshotCenter(s,visibilities,times)
    MapMats.inverseCovarianceWeightVisibilities(s,visibilities)

    
    #Perform mapmaking and calculate PSFs
    coaddedMap = np.zeros(coords.nFacetPixels)
    PSF = np.zeros((coords.nFacetPixels,coords.nExtendedPixels))
    pointSourcePSF = np.zeros((coords.nFacetPixels, ps.nSources))
    for snapshot in times.snapshots:    
        NinvTimesy = MapMats.calculateNinvTimesy(visibilities, snapshot)
        Ninv = MapMats.calculateNInv(s,snapshot)
        KAtranspose = MapMats.calculateKAtranspose(s,snapshot,coords,PBs)    
        coaddedMap += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], NinvTimesy))
        PSF += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], np.dot(np.diag(Ninv), KAtranspose.conj().transpose())))    
        if s.PSFforPointSources and ps.nSources > 0:
            pointSourceAmatrix = MapMats.calculatePSAmatrix(s,snapshot,ps,PBs)
            pointSourcePSF += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], np.dot(np.diag(Ninv), pointSourceAmatrix)))
            
    #Renormalize maps and PSFs and save results
    Dmatrix = np.diag(np.diag(PSF[:,coords.mapIndexLocationsInExtendedIndexList])**(-1))
    #Dmatrix = np.diag(np.ones(len(coords.mapIndices)))
    #print "warning: Dmatrix set to Identity"
    PSF = np.dot(Dmatrix,PSF)
    coaddedMap = np.dot(Dmatrix,coaddedMap)
    mapNoiseCovariance = np.dot(PSF[:,coords.mapIndexLocationsInExtendedIndexList],np.transpose(Dmatrix))
    if s.PSFforPointSources and ps.nSources > 0:
        pointSourcePSF = np.dot(Dmatrix,pointSourcePSF)
    
    


    s.GSMNSIDE = s.mapNSIDE
    coordsGSM = Geometry.Coordinates(s,True)
    GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
    interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
    interpoltedGSMRotated = np.zeros(len(interpoltedGSMRotated))
    testSourceIndex = hp.query_disc(coords.NSIDE, hp.ang2vec(np.pi/2, 0), s.facetSize/10 * 2*np.pi/360.0)[0]
    interpoltedGSMRotated[testSourceIndex] = 1
    convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])
    plt.figure()    
    hp.mollview(interpoltedGSMRotated, title="interpoltedGSMRotated")
    
    def plotFacet(s,coords,facetMap,plotTitle):
        mapToPlot = np.zeros(coords.mapPixels)
        mapToPlot[coords.mapIndices] = facetMap    
        hp.mollview(mapToPlot, title=plotTitle)
        plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)
    plt.figure()
    plt.plot(visibilities)
    plotFacet(s,coords,interpoltedGSMRotated[coords.mapIndices],"GSM")    

    mapToPlot = np.zeros(coords.mapPixels)
    mapToPlot[coords.mapIndices] = PSF[:,np.nonzero(coords.extendedIndices == testSourceIndex)[0][0]]    
    hp.mollview(mapToPlot, title="PSF")
    plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)


    plotFacet(s,coords,convolvedGSM,"Convolved GSM")
    plotFacet(s,coords,coaddedMap,"coaddedMap")
    plotFacet(s,coords,coaddedMap/convolvedGSM,'map/convolved ratio')
      
    plt.plot(np.diag(Dmatrix))
    plt.show()


    
    
    
    
    
    MapMats.saveAllResults(s,coords,times,ps,Dmatrix,PSF,coaddedMap,pointSourcePSF,mapNoiseCovariance)
    return s.resultsFolder
    
    

if __name__ == "__main__":
    Mapmaker()