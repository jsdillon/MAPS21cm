import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
import ephem
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
from Source.GlobalSkyModel import GlobalSkyModel
import scipy.constants as const
import cPickle as pickle
import os
from Mapmaker import Mapmaker

plt.close("all")

def loadAllResults(resultsFolder):
    s = pickle.load(open(resultsFolder + "specifications.p","rb"))
    coords = pickle.load(open(s.resultsFolder + "coordinates.p","rb"))
    times = pickle.load( open(s.resultsFolder + "times.p","rb"))
    ps = pickle.load(open(s.resultsFolder + "pointSourceCatalog.p","rb"))
    Dmatrix = np.load(s.resultsFolder + "Dmatrix.npy")
    PSF = np.load(s.resultsFolder + "PSF.npy")
    coaddedMap = np.load(s.resultsFolder + "coaddedMap.npy")
    mapNoiseCovariance = np.load(s.resultsFolder + "mapNoiseCovariance.npy")
    if s.PSFforPointSources and ps.nSources>0:
        pointSourcePSF = np.load(s.resultsFolder + "pointSourcePSF.npy")
        return s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF
    else:
        return s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, []
    
def plotFacet(s,coords,facetMap,plotTitle):
    mapToPlot = np.zeros(coords.mapPixels)
    mapToPlot[coords.mapIndices] = facetMap    
    hp.mollview(mapToPlot, title=plotTitle)
    plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)

    
#Test 1: GSM Only    
def TestGSMOnly():
    print "\nNow running PSF/Mapmaking Comparison GSM Only..."
    resultsDirectory = Mapmaker(PSFextensionBeyondFacetFactor = 2, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = False)
    s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = loadAllResults(resultsDirectory)
    s.GSMNSIDE = s.mapNSIDE
    coordsGSM = Geometry.Coordinates(s,True)
    GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
    interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
    convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])
    #plotFacet(s,coords,convolvedGSM,"Convolved GSM")
    #plotFacet(s,coords,coaddedMap,"Coadded Map")
    #plotFacet(s,coords,coaddedMap/convolvedGSM,"Coadded Map / Convolved GSM")
    print "Error = " + str(np.linalg.norm(coaddedMap - convolvedGSM)/np.linalg.norm(convolvedGSM))

#Test 2: Point Sources Only
def TestPointSourcesOnly():
    print "\nNow running PSF/Mapmaking Comparison Point Sources Only..."
    resultsDirectory = Mapmaker(simulateVisibilitiesWithGSM = False, simulateVisibilitiesWithPointSources = True)
    s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = loadAllResults(resultsDirectory)
    convolvedPointSources = np.dot(pointSourcePSF, ps.scaledFluxes)
    #plotFacet(s,coords,convolvedPointSources,"Convolved Sources")
    #plotFacet(s,coords,coaddedMap,"Coadded Map")
    #plotFacet(s,coords,coaddedMap/convolvedPointSources,"Coadded Map / Convolved Point Sources")
    print "Error = " + str(np.linalg.norm(coaddedMap - convolvedPointSources)/np.linalg.norm(convolvedPointSources))

#Test 3: Error as a function of PSFextensionBeyondFacetFactor
def TestErrorVsPSFext():
    PSFErrors = []
    PSFextFactors = np.arange(1,2.6,.1)
    for extFactor in PSFextFactors:
        resultsDirectory = Mapmaker(PSFextensionBeyondFacetFactor = extFactor, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True)
        s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = loadAllResults(resultsDirectory)
        convolvedPointSources = np.dot(pointSourcePSF, ps.scaledFluxes)
        s.GSMNSIDE = s.mapNSIDE
        coordsGSM = Geometry.Coordinates(s,True)
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
        convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])
        PSFErrors.append(np.linalg.norm(coaddedMap - convolvedPointSources - convolvedGSM)/np.linalg.norm(convolvedPointSources + convolvedGSM))
    plt.semilogy(PSFextFactors, PSFErrors)
    plt.xlabel('PSF Size Relative to Facet Size')
    plt.ylabel('PSF Convolution vs. Mapmaking Error')

#Test 4: Error as a function of integrationsPerSnapshot
def TestErrorVsIntegrations():
    resultsDirectory = Mapmaker(PSFextensionBeyondFacetFactor = 3, integrationsPerSnapshot = 1, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True, MaximumAllowedAngleFromFacetCenterToPointingCenter = 5, GSMNSIDE = 32, mapNSIDE = 32)
    s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = loadAllResults(resultsDirectory)
    trueMap = np.copy(coaddedMap)
    intErrors = []
    integrationsList = np.arange(3,30,2)
    for integrations in integrationsList:
        resultsDirectory = Mapmaker(PSFextensionBeyondFacetFactor = 3, integrationsPerSnapshot = integrations, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True, MaximumAllowedAngleFromFacetCenterToPointingCenter = 5, GSMNSIDE = 32, mapNSIDE = 32)
        s, coords, times, ps, Dmatrix, PSF, coaddedMap, mapNoiseCovariance, pointSourcePSF = loadAllResults(resultsDirectory)
        intErrors.append(np.linalg.norm(coaddedMap - trueMap)/np.linalg.norm(trueMap))
    plt.semilogy(integrationsList, intErrors)
    plt.xlabel('Number of Integrations Per Snapshot')
    plt.ylabel('Error Compared to 1 Integration Per Snapshot')


#TestGSMOnly()
#TestPointSourcesOnly()
#TestErrorVsPSFext()
TestErrorVsIntegrations()
