# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ephem
import Geometry
from GlobalSkyModel import GlobalSkyModel


def VisibilitySimulator(s,PBs,ps,times):
    print "Now simulating visibilities (assuming XX beams only)..."
    coordsGSM = Geometry.Coordinates(s,True)
    visibilities = np.zeros([len(times.LSTs),len(s.baselines)],dtype=complex)
    
    #TODO: this ignores polarization and differing primary beams
    if s.simulateVisibilitiesWithGSM:
        GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
        #interpolate onto rotated equatorial coordinates
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
        hp.mollview(np.log10(interpoltedGSMRotated), title="GSM in Rotated Equatorial Coordinates")
        
        #loop over times and baselines to calculate visibilities
        for t in range(len(times.LSTs)):
            pixelAlts, pixelAzs = Geometry.convertEquatorialToHorizontal(s,coordsGSM.pixelRAs,coordsGSM.pixelDecs,times.LSTs[t])
            rHatVectors = Geometry.convertAltAzToCartesian(pixelAlts,pixelAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[t]), np.pi/2-pixelAlts, pixelAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-2j*np.pi*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(GSM.hpMap * primaryBeam * exponent) * 4*np.pi / len(GSM.hpMap) / s.convertJyToKFactor

    if s.simulateVisibilitiesWithPointSources:
        for t in range(len(times.LSTs)):
            psAlts, psAzs = Geometry.convertEquatorialToHorizontal(s,ps.RAs,ps.decs,times.LSTs[t])
            rHatVectors = Geometry.convertAltAzToCartesian(psAlts,psAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[t]), np.pi/2-psAlts, psAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-2j*np.pi*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(ps.scaledFluxes * primaryBeam * exponent) 
				
    return visibilities