# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13 as cosmo

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
        
        if s.SquareFacetInRADec:
            deltaDecs = np.abs(self.pixelDecs - s.facetDecinRad)     
            deltaRAs = np.minimum(np.minimum(np.abs(self.pixelRAs - s.facetRAinRad), np.abs(self.pixelRAs - s.facetRAinRad + 2*np.pi)), np.abs(self.pixelRAs - s.facetRAinRad - 2*np.pi))
            self.mapIndices = np.flatnonzero((deltaDecs <= s.facetSize * 2*np.pi/360.0/2.0) * (deltaRAs <= s.facetSize * 2*np.pi/360.0/2.0/np.cos(s.facetDecinRad)))
            self.extendedIndices = np.nonzero((deltaDecs <= s.PSFextensionBeyondFacetFactor*s.facetSize * 2*np.pi/360.0/2.0) * (deltaRAs <= s.PSFextensionBeyondFacetFactor*s.facetSize * 2*np.pi/360.0/2.0/np.cos(s.facetDecinRad)))[0]
        else:
            self.mapIndices = hp.query_disc(self.NSIDE, hp.ang2vec(np.pi/2, 0), s.facetSize * 2*np.pi/360.0)
            self.extendedIndices = hp.query_disc(self.NSIDE, hp.ang2vec(np.pi/2, 0), s.facetSize * 2*np.pi/360.0 * s.PSFextensionBeyondFacetFactor)        
        self.nFacetPixels = len(self.mapIndices)
        self.nExtendedPixels = len(self.extendedIndices)
        
        extendedIndexDict = dict([ (self.extendedIndices[i], i) for i in range(self.nExtendedPixels) ])
        self.mapIndexLocationsInExtendedIndexList = np.asarray([extendedIndexDict[mapIndex] for mapIndex in self.mapIndices])
        self.mapIndexOfFacetCenter = np.dot(self.mapIndices ==  hp.ang2pix(self.NSIDE,np.pi/2,0.0),np.arange(len(self.mapIndices)))
        self.extendedIndexOfFacetCenter = np.dot(self.extendedIndices ==  hp.ang2pix(self.NSIDE,np.pi/2,0.0),np.arange(len(self.extendedIndices)))
        
    def computeCubeCoordinates(self,s,freqs):
        """ Given a range of frequencies, this function calculates the coordinates of each voxel in cMpc as a function of [freq,pixelIndex] and other useful quantities."""
        self.freqs = freqs       
        self.nFreqs = len(freqs)
        self.nVoxels = self.nFacetPixels * self.nFreqs
        self.comovingDistances = cosmo.comoving_distance(1420.40575177 / freqs - 1).value
        self.comovingTransverseDistances = cosmo.comoving_distance(1420.40575177 / freqs - 1).value
        self.losCoords = np.outer(self.comovingDistances,np.ones(self.nFacetPixels)) #distance from the earth
        self.xCoords = np.outer(self.comovingTransverseDistances,hp.pix2vec(self.NSIDE, self.mapIndices)[1]) #distance from the facet center, roughly the same direction as RA
        self.yCoords = np.outer(self.comovingTransverseDistances,hp.pix2vec(self.NSIDE, self.mapIndices)[2]) #distance from the facet center, roughly the same direction as Dec
        self.losDelta = np.abs(np.mean(np.diff(self.comovingDistances)))
        self.cubeCoordinates = np.zeros((len(freqs),self.nFacetPixels), dtype=('f8,f8,f8'))


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


class Snapshot:
    """ This class has information about a single snapshot. """
    def __init__(self, times, LSTindices):
        self.LSTindices = LSTindices
        self.LSTs = times.LSTs[LSTindices]
        self.centralLSTIndex = LSTindices[int(round(len(self.LSTindices)/2.0 - .5))]
        self.centralLST = times.LSTs[self.centralLSTIndex]

class Times:
    """ This class contains information about the integration and snapshot LSTs. """
    def __init__(self, s):   
        self.LSTs = np.loadtxt(s.LSTsFilename)
        self.integrationTime = np.median(self.LSTs[1:] - self.LSTs[0:len(self.LSTs)-1]) * 60 * 60
        LSTrange = self.LSTs[-1]*360.0/24 - self.LSTs[0]*360.0/24
        if np.abs(s.facetRA - self.LSTs[0]*360.0/24) < LSTrange/4 or np.abs(s.facetRA - self.LSTs[-1]*360.0/24) < LSTrange/4: #this fixes the problem when the facet center is near RA=0
            self.LSTs = np.fft.fftshift(self.LSTs)
            s.pointings = np.fft.fftshift(s.pointings, axes=0) #only shifts along the LST axis, not the antenna axis
            s.noisePerAntenna = np.fft.fftshift(s.noisePerAntenna, axes=0) #only shifts along the LST axis, not the antenna axis
        self.useThisLST = np.ones(len(self.LSTs))
        self.snapshots = [] 
        
    # Calculates which LSTs are close enough to the facet center to be used in making the map, eliminating the rest from loaded files 
    def CutOutUnusedLSTsAndGroupIntoSnapshots(self, s):
        for t in range(len(self.LSTs)):
            facetCenterAlt, facetCenterAz = convertEquatorialToHorizontal(s, s.facetRAinRad, s.facetDecinRad, self.LSTs[t])        
            if s.antennasHaveIdenticalBeams:
                angularDistanceFromFacetToPointing = hp.rotator.angdist([np.pi/2 - facetCenterAlt, facetCenterAz],[np.pi/2 - s.pointingCenters[s.pointings[t]][0], s.pointingCenters[s.pointings[t]][1]])
            else:
                angularDistanceFromFacetToPointing = np.min(np.array([hp.rotator.angdist([np.pi/2 - facetCenterAlt, facetCenterAz],[np.pi/2 - s.pointingCenters[s.pointings[t,ant]][0], s.pointingCenters[s.pointings[t,ant]][1]]) for ant in range(s.nAntennas)]))
            if angularDistanceFromFacetToPointing > s.MaximumAllowedAngleFromFacetCenterToPointingCenter*2.0*np.pi/360.0:
                self.useThisLST[t] = False

        #Groups LST indices into snapshots, assuming that all LSTs are evenly spaced in time
        LSTindicesToUse = np.nonzero(self.useThisLST)[0]
        closeEnoughLSTs = len(LSTindicesToUse)
        consecutiveSegments =  np.split(LSTindicesToUse, np.where(np.diff(LSTindicesToUse) != 1)[0]+1)
        LSTindicesToUse = np.asarray([segment[index] for segment in consecutiveSegments for index in range(int(np.floor(len(segment)/s.integrationsPerSnapshot)*s.integrationsPerSnapshot))])
        
        #Delete unnecessary data
        self.useThisLST = np.zeros(len(self.LSTs))
        self.useThisLST[LSTindicesToUse] = True
        self.LSTs = self.LSTs[self.useThisLST == True]
        s.pointings = s.pointings[self.useThisLST == True]
        s.noisePerAntenna = s.noisePerAntenna[self.useThisLST == True]
            
        #Make snapshots
        LSTindicesGrouped = np.reshape(np.arange(len(self.LSTs)),(len(LSTindicesToUse)/s.integrationsPerSnapshot, s.integrationsPerSnapshot))
        self.snapshots = [Snapshot(self, LSTgroup) for LSTgroup in LSTindicesGrouped]        
        print "Observations of " + str(closeEnoughLSTs) + " LSTs are within " + str(s.MaximumAllowedAngleFromFacetCenterToPointingCenter) + " degrees of the facet center."
        print "They are broken up into " + str(len(self.snapshots)) + " snapshots of exactly " + str(s.integrationsPerSnapshot * self.integrationTime) + " seconds.\n" + str(closeEnoughLSTs - len(self.LSTs)) + " integration(s) were discarded in snapshotting."
    
def rephaseVisibilitiesToSnapshotCenter(s,visibilities,times):
    """This function steps through the snapshots and rephases all the visibilities to the facet center of central LST of the snapshot."""
    for snapshot in times.snapshots:
        centralLSTAlt, centralLSTaz = convertEquatorialToHorizontal(s,s.facetRAinRad,s.facetDecinRad,snapshot.centralLST)
        centralSnapshotFacetCenterVector = convertAltAzToCartesian(centralLSTAlt, centralLSTaz)        
        for t in range(len(snapshot.LSTs)):
            thisLSTAlt, thisLSTaz = convertEquatorialToHorizontal(s,s.facetRAinRad,s.facetDecinRad,snapshot.LSTs[t])
            deltaTheta = centralSnapshotFacetCenterVector - convertAltAzToCartesian(thisLSTAlt, thisLSTaz)
            rephaseFactors = np.exp(-1j * s.k * s.baselines.dot(deltaTheta))
            visibilities[snapshot.LSTindices[t],:] = visibilities[snapshot.LSTindices[t],:] * rephaseFactors

