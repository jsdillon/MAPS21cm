import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from MAPS21cm.Specifications import Specifications
from MAPS21cm.PrimaryBeams import PrimaryBeams
from MAPS21cm.VisibilitySimulator import VisibilitySimulator
import ephem
from MAPS21cm import Geometry
from MAPS21cm.PointSourceCatalog import PointSourceCatalog
from MAPS21cm import MatricesForMapmaking as MapMats
from MAPS21cm.GlobalSkyModel import GlobalSkyModel
import scipy.constants as const
import cPickle as pickle
import os
from MAPS21cm.Mapmaker import Mapmaker

plt.close("all")

def plotFacet(s,coords,facetMap,plotTitle):
    mapToPlot = np.zeros(coords.mapPixels)
    mapToPlot[coords.mapIndices] = facetMap    
    hp.mollview(mapToPlot, title=plotTitle)
    plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)

mainDirectory = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0] #Directory above this one

    
#Test 1: GSM Only    
def TestGSMOnly():
    print "\nNow running PSF/Mapmaking Comparison GSM Only..."
    resultsDirectory = Mapmaker(mainDirectory, PSFextensionBeyondFacetFactor = 8, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = False)
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    if s.GSMNSIDE < s.mapNSIDE:
        s.GSMNSIDE = s.mapNSIDE
    if s.GSMNSIDE > s.mapNSIDE:
        "WARNING: This test will fail because GSMNSIDE > mapNSIDE."
    coordsGSM = Geometry.Coordinates(s,True)
    coords = Geometry.Coordinates(s)
    print "There are " + str(coords.nExtendedPixels) + " pixels in the extended facet."
    GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
    interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
    convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])
    plotFacet(s,coords,convolvedGSM,"Convolved GSM")
    plotFacet(s,coords,coaddedMap,"Coadded Map")
    plotFacet(s,coords,coaddedMap/convolvedGSM,"Coadded Map / Convolved GSM")
    print "Error = " + str(np.linalg.norm(coaddedMap - convolvedGSM)/np.linalg.norm(convolvedGSM))

#Test 2: Point Sources Only
def TestPointSourcesOnly():
    print "\nNow running PSF/Mapmaking Comparison Point Sources Only..."
    resultsDirectory = Mapmaker(mainDirectory, simulateVisibilitiesWithGSM = False, simulateVisibilitiesWithPointSources = True)
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    convolvedPointSources = np.dot(pointSourcePSF, ps.scaledFluxes)
    #plotFacet(s,coords,convolvedPointSources,"Convolved Sources")
    #plotFacet(s,coords,coaddedMap,"Coadded Map")
    #plotFacet(s,coords,coaddedMap/convolvedPointSources,"Coadded Map / Convolved Point Sources")
    print "Error = " + str(np.linalg.norm(coaddedMap - convolvedPointSources)/np.linalg.norm(convolvedPointSources))

#Test 3: Error as a function of PSFextensionBeyondFacetFactor
def TestErrorVsPSFext():    
    PSFErrors = []
    PSFextFactors = np.arange(1,10,.5)
    for extFactor in PSFextFactors:
        resultsDirectory = Mapmaker(mainDirectory, PSFextensionBeyondFacetFactor = extFactor, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True)
        s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
        convolvedPointSources = np.dot(pointSourcePSF, ps.scaledFluxes)
        if s.GSMNSIDE < s.mapNSIDE:
            s.GSMNSIDE = s.mapNSIDE
        coordsGSM = Geometry.Coordinates(s,True)
        coords = Geometry.Coordinates(s)        
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
        convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])
        PSFErrors.append(np.linalg.norm(coaddedMap - convolvedPointSources - convolvedGSM)/np.linalg.norm(convolvedPointSources + convolvedGSM))
    plt.figure()
    plt.semilogy(PSFextFactors, PSFErrors)
    plt.xlabel('PSF Size Relative to Facet Size')
    plt.ylabel('PSF Convolution vs. Mapmaking Error')

#Test 4: Error as a function of integrationsPerSnapshot
def TestErrorVsIntegrations():
    resultsDirectory = Mapmaker(mainDirectory, PSFextensionBeyondFacetFactor = 3, integrationsPerSnapshot = 1, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True, MaximumAllowedAngleFromFacetCenterToPointingCenter = 4.33, GSMNSIDE = 32, mapNSIDE = 32, facetRA = .01)
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    coords = Geometry.Coordinates(s)    
    #plotFacet(s,coords,coaddedMap,'Coadded Map with 1 Integration per Snapshot')    
    trueMap = np.copy(coaddedMap)
    intErrors = []
    integrationsList = np.asarray([3, 9, 27, 81, 243])
    for integrations in integrationsList:
        resultsDirectory = Mapmaker(mainDirectory, PSFextensionBeyondFacetFactor = 3, integrationsPerSnapshot = integrations, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = True, MaximumAllowedAngleFromFacetCenterToPointingCenter = 4.33, GSMNSIDE = 32, mapNSIDE = 32, facetRA = .01)
        s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
        coords = Geometry.Coordinates(s)        
        #plotFacet(s,coords,coaddedMap,'Coadded Map with ' + str(integrations) + ' Integrations per Snapshot')
        intErrors.append(np.linalg.norm(coaddedMap - trueMap)/np.linalg.norm(trueMap))
    plt.figure()    
    plt.loglog(integrationsList, intErrors)
    plt.xlabel('Number of Integrations Per Snapshot')
    plt.ylabel('Error Compared to 1 Integration Per Snapshot')


###############################################################################################################################
#   VARIOUS TESTS OF THE MAPMAKING ALGORITHM REPRODUCING DILLON ET AL. (2015) RESULTS
###############################################################################################################################

TestGSMOnly()
TestPointSourcesOnly()
TestErrorVsPSFext()
TestErrorVsIntegrations()

