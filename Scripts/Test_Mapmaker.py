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
from astropy import units as u
from astropy.coordinates import SkyCoord

plt.close("all")

def plotFacet(s,coords,facetMap,plotTitle):
    if s.makeFacetSameAsAdaptivePSF and s.useAdaptiveHEALPixForPSF:
        hp.mollview(facetMap[coords.newPSFIndices], title=plotTitle)
    else:
        mapToPlot = np.zeros(coords.mapPixels)
        mapToPlot[coords.facetIndices] = facetMap    
        hp.mollview(mapToPlot, title=plotTitle)
        plt.axis(np.asarray([-1,1,-1,1]) * s.facetSize/50)

mainDirectory = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0] #Directory above this one

def AdaptiveResolutionGSM(s, coords):
    """This function computes the GSM for the adaptive PSF where each pixel comes from a map at the appropriate resolution"""
    galCoords = SkyCoord(frame="icrs", ra=coords.PSFRAs*u.rad, dec=coords.PSFDecs*u.rad).transform_to("galactic")
    NSIDE = s.adaptiveHEALPixMinNSIDE
    interpoltedGSMRotated = np.zeros(coords.nPSFPixels)    
    while NSIDE <= s.mapNSIDE:
        thisGSM = GlobalSkyModel(s.freq, s.GSMlocation, NSIDE).hpMap
        interpoltedGSMRotated[coords.newPSFNSIDEs==NSIDE] = hp.get_interp_val(thisGSM,-galCoords[coords.newPSFNSIDEs==NSIDE].b.radian+np.pi/2, galCoords[coords.newPSFNSIDEs==NSIDE].l.radian)
        NSIDE *= 2

    return interpoltedGSMRotated
    
#Test 1: GSM Only    
def TestGSMOnly():
    print "\nNow running PSF/Mapmaking Comparison GSM Only..."
    #resultsDirectory = Mapmaker(mainDirectory, PSFextensionBeyondFacetFactor = 8, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = False)
    resultsDirectory = Mapmaker(mainDirectory, simulateVisibilitiesWithGSM = True, simulateVisibilitiesWithPointSources = False)
    s, times, ps, Dmatrix, PSF, coaddedMap, pointSourcePSF = MapMats.loadAllResults(resultsDirectory)
    if s.GSMNSIDE < s.mapNSIDE:
        s.GSMNSIDE = s.mapNSIDE
    if s.GSMNSIDE > s.mapNSIDE:
        "WARNING: This test will fail because GSMNSIDE > mapNSIDE."
    coordsGSM = Geometry.Coordinates(s,useAnotherResolution = s.GSMNSIDE)
    coords = Geometry.Coordinates(s)
    if s.useAdaptiveHEALPixForPSF: coords.convertToAdaptiveHEALPix(s, times)
    print "There are " + str(coords.nPSFPixels) + " pixels in the PSF."
      
    if s.useAdaptiveHEALPixForPSF:    
        interpoltedGSMRotated = AdaptiveResolutionGSM(s, coords)
    else:
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        galCoords = SkyCoord(frame="icrs", ra=coords.PSFRAs*u.rad, dec=coords.PSFDecs*u.rad).transform_to("galactic")        
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-galCoords.b.radian+np.pi/2, galCoords.l.radian)
    convolvedGSM = np.dot(PSF,interpoltedGSMRotated)
    
    plotFacet(s,coords,convolvedGSM,"Convolved GSM")
    plotFacet(s,coords,coaddedMap,"Coadded Map")
    plotFacet(s,coords,coaddedMap - convolvedGSM,"Coadded Map - Convolved GSM")
    print "Error = " + str(np.linalg.norm(coaddedMap - convolvedGSM)/np.linalg.norm(convolvedGSM))
    
#    plt.figure()
#    plt.plot(coords.PSFDecs,'.')
#    plt.figure()
#    #plt.scatter(coords.PSFRAs, coords.PSFDecs, c = np.log10(np.abs(PSF[coords.facetIndexOfFacetCenter,:])), s = 1024 / coords.newPSFNSIDEs, edgecolor='none')
#    plt.scatter(coords.PSFRAs, coords.PSFDecs, c = PSF[coords.facetIndexOfFacetCenter,:], s = 1024 / coords.newPSFNSIDEs, edgecolor='none')
#    plt.colorbar()
#    plt.clim(-6, 0)

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
        coordsGSM = Geometry.Coordinates(s,useAnotherResolution = s.GSMNSIDE)
        coords = Geometry.Coordinates(s)        
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
        convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.PSFIndices])
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
#TestPointSourcesOnly()
#TestErrorVsPSFext()
#TestErrorVsIntegrations()

