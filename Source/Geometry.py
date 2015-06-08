# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt


class Coordinates:
    """This class figures out the RA, Dec, and Galactic Coordinates of every pixel in the map where the facet center is rotated to lie on the horizon. Detaults to the map resolution, but can also be done for the GSM resolution used."""
    def __init__(self, s, atResolutionForSimulation = False):
        if atResolutionForSimulation:
            self.NSIDE = s.GSMNSIDE
        else:
            self.NSIDE = s.mapNSIDE
        self.mapPixels = 12 * self.NSIDE**2

        healpixCoords = np.transpose(np.asarray(hp.pix2ang(self.NSIDE,np.arange(self.mapPixels))))                
        self.originalPixelRAs = healpixCoords[:,1]
        self.originalPixelDecs = -healpixCoords[:,0] + np.pi/2
        self.originalGalCoords = SkyCoord(frame="icrs", ra=self.originalPixelRAs*u.rad, dec=self.originalPixelDecs*u.rad).transform_to("galactic")                
        originalPixelVectors = [ np.array([np.cos(self.originalPixelDecs[n]) * np.cos(self.originalPixelRAs[n]), np.cos(self.originalPixelDecs[n]) * np.sin(self.originalPixelRAs[n]), np.sin(self.originalPixelDecs[n])]) for n in range(self.mapPixels)]
        
        rotateToPrimeMeridian = np.array([[np.cos(s.facetRAinRad), -np.sin(s.facetRAinRad), 0], [np.sin(s.facetRAinRad), np.cos(s.facetRAinRad), 0],[0, 0, 1]])
        rotateToEquator = np.array([[np.cos(-s.facetDecinRad), 0, np.sin(-s.facetDecinRad)], [0, 1, 0], [-np.sin(-s.facetDecinRad), 0, np.cos(-s.facetDecinRad)]])
        rotatedVectors = [np.dot(rotateToPrimeMeridian, np.dot(rotateToEquator,originalPixelVectors[n])) for n in range(self.mapPixels)]
        
        self.pixelDecs = np.asarray([np.arcsin(rotatedVectors[n][2]) for n in range(self.mapPixels)])
        self.pixelRAs = np.asarray([np.arctan2(rotatedVectors[n][1],rotatedVectors[n][0]) for n in range(self.mapPixels)])
        self.pixelRAs[self.pixelRAs < 0] = self.pixelRAs[self.pixelRAs < 0] + 2*np.pi
        self.galCoords = SkyCoord(frame="icrs", ra=self.pixelRAs*u.rad, dec=self.pixelDecs*u.rad).transform_to("galactic")

                        
# Convert RAs and Decs in radians to altitudes and azimuths in radians, given an LST and array location in the specs object
def convertEquatorialToHorizontal(s,RAs,decs,LST):
    """ Convert RAs and Decs in radians to altitudes and azimuths in radians, given an LST and array location in the specs object """
    lha = np.pi/12.0*LST - RAs #in radians
    altitudes = np.arcsin(np.sin(s.arrayLatInRad) * np.sin(decs) + np.cos(s.arrayLatInRad) * np.cos(decs) * np.cos(lha))
    azimuths = np.arctan2( np.sin(lha) * np.cos(decs), np.cos(lha) * np.cos(decs) * np.sin(s.arrayLatInRad) - np.sin(decs) * np.cos(s.arrayLatInRad)) + np.pi;	
    return altitudes, azimuths


def convertAltAzToCartesian(alts, azs):
    """ Convert list of altitudes and azimuths to cartesian coordinates."""        
    return np.transpose(np.asarray([np.sin(azs)*np.cos(alts), np.cos(azs)*np.cos(alts), np.sin(alts)]))

# Calculates which LSTs are close enough to the facet center to be used in making the map, eliminating the rest from loaded files 
def CutOutUnusedLSTs(s):
    for t in range(len(s.LSTs)):
        facetCenterAlt, facetCenterAz = convertEquatorialToHorizontal(s, s.facetRAinRad, s.facetDecinRad, s.LSTs[t])        
        if s.antennasHaveIdenticalBeams:
            angularDistanceFromFacetToPointing = hp.rotator.angdist([np.pi/2 - facetCenterAlt, facetCenterAz],[np.pi/2 - s.pointingCenters[s.pointings[t]][0], s.pointingCenters[s.pointings[t]][1]])
        else:
            angularDistanceFromFacetToPointing = np.min(np.array([hp.rotator.angdist([np.pi/2 - facetCenterAlt, facetCenterAz],[np.pi/2 - s.pointingCenters[s.pointings[t,ant]][0], s.pointingCenters[s.pointings[t,ant]][1]]) for ant in range(s.nAntennas)]))
        if angularDistanceFromFacetToPointing > s.MaximumAllowedAngleFromFacetCenterToZenith*2.0*np.pi/360.0:
            s.useThisLST[t] = False

    s.LSTs = s.LSTs[s.useThisLST == True]
    s.pointings = s.pointings[s.useThisLST == True]
    if s.useOnlyUniqueBaselines:
        s.noisePerUniqueBaseline = s.noisePerUniqueBaseline[s.useThisLST == True]
    else:
        s.noisePerAntenna = s.noisePerAntenna[s.useThisLST == True]
    print "Observations of " + str(int(np.sum(s.useThisLST))) + " LSTs are within " + str(s.MaximumAllowedAngleFromFacetCenterToZenith) + " degrees of the facet center."


