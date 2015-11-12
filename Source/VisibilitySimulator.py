# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ephem
import Geometry
from GlobalSkyModel import GlobalSkyModel

def AdaptiveResolutionGSM(s, coords):
    """This function computes the GSM for the adaptive PSF where each pixel comes from a map at the appropriate resolution"""
    from astropy import units as u
    from astropy.coordinates import SkyCoord    
    galCoords = SkyCoord(frame="icrs", ra=coords.PSFRAs*u.rad, dec=coords.PSFDecs*u.rad).transform_to("galactic")
    NSIDE = s.adaptiveHEALPixMinNSIDE
    interpoltedGSMRotated = np.zeros(coords.nPSFPixels)    
    while NSIDE <= s.mapNSIDE:
        thisGSM = GlobalSkyModel(s.freq, s.GSMlocation, NSIDE).hpMap
        interpoltedGSMRotated[coords.newPSFNSIDEs==NSIDE] = hp.get_interp_val(thisGSM,-galCoords[coords.newPSFNSIDEs==NSIDE].b.radian+np.pi/2, galCoords[coords.newPSFNSIDEs==NSIDE].l.radian)
        NSIDE *= 2
    return interpoltedGSMRotated

def VisibilitySimulator(s,PBs,ps,times,coords):
    print "Now simulating visibilities (assuming XX beams only)..."
    if s.GSMNSIDE < s.mapNSIDE:
        s.GSMNSIDE = s.mapNSIDE
    coordsGSM = Geometry.Coordinates(s,useAnotherResolution = s.GSMNSIDE)
    visibilities = np.zeros([len(times.LSTs),len(s.baselines)],dtype=complex)
    
    #TODO: this ignores polarization and differing primary beams
    if s.simulateVisibilitiesWithGSM:
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        #print "using adaptive GSM in visibility simulation"
        #interpoltedGSMRotated = AdaptiveResolutionGSM(s, coords)
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))        
        
        #loop over times and baselines to calculate visibilities
        for t in range(len(times.LSTs)):
            pixelAlts, pixelAzs = Geometry.convertEquatorialToHorizontal(s,coordsGSM.pixelRAs,coordsGSM.pixelDecs,times.LSTs[t])
            #pixelAlts, pixelAzs = Geometry.convertEquatorialToHorizontal(s,coords.PSFRAs,coords.PSFDecs,times.LSTs[t])
            rHatVectors = Geometry.convertAltAzToCartesian(pixelAlts,pixelAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[t]), np.pi/2-pixelAlts, pixelAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-1j*s.k*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(interpoltedGSMRotated * primaryBeam * exponent) * 4*np.pi / len(GSM.hpMap) / s.convertJyToKFactor
                #visibilities[t,b] += np.sum(interpoltedGSMRotated * primaryBeam * exponent * 4*np.pi / (12 * coords.newPSFNSIDEs**2)) / s.convertJyToKFactor

    if s.simulateVisibilitiesWithPointSources and ps.nSources > 0:
        for t in range(len(times.LSTs)):
            psAlts, psAzs = Geometry.convertEquatorialToHorizontal(s,ps.RAs,ps.decs,times.LSTs[t])
            rHatVectors = Geometry.convertAltAzToCartesian(psAlts,psAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[t]), np.pi/2-psAlts, psAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-1j*s.k*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(ps.scaledFluxes * primaryBeam * exponent) 
				
    return visibilities