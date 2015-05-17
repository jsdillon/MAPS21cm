# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ephem
import Source.MapmakerHelper as mmh
from Source.GlobalSkyModel import GlobalSkyModel as gsm


def VisibilitySimulator(s,g,PBs):
    print "Now simulating visibilities..."
    visibilities = np.zeros([len(s.LSTs),len(s.baselines)],dtype=complex)
    
    #TODO: this ignores polarization and differing primary beams
    if s.simulateVisibilitiesWithGSM:
        GSM = gsm(s.freq, s.GSMlocation, s.GSMNSIDE)
        #interpolate onto rotated equatorial coordinates
        interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-g.galCoords.b.radian+np.pi/2, np.asarray(g.galCoords.l.radian))
        hp.mollview(np.log10(interpoltedGSMRotated), title="GSM in Rotated Equatorial Coordinates")


#TODO: investigate temperature/flux normalization        
        
        
#TODO: FIX WHAT HAPPENS WHEN GSM AND MAP RESOLUTION AREN'T THE SAME RESOLUTION  
        
        for t in range(len(s.LSTs)):
            pixelAlts, pixelAzs = mmh.convertEquatorialToHorizontal(s,g.pixelRAs,g.pixelDecs,s.LSTs[t])
            rHatVectors = mmh.convertAltAzToCartesian(pixelAlts,pixelAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[0]), np.pi/2-pixelAlts, pixelAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-2j*np.pi*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(GSM.hpMap * primaryBeam * exponent) * 4*np.pi / len(GSM.hpMap)
                
#    if s.simulateVisibilitiesWithPointSources:
        #Load PS catalog
        

    return visibilities