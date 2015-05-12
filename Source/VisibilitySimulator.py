# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ephem
import Source.MapmakerHelper as mmh


def VisibilitySimulator(s,g,PBs):
    print "Now simulating visibilities..."
    visibilities = np.zeros([len(s.LSTs),len(s.baselines)],dtype=complex)
    
    #TODO: this ignores polarization and differing primary beams
    if s.simulateVisibilitiesWithGSM:
        #Compute GSM from 3 principal components appropriately weighted and then interpolate onto rotated equatorial coordinates
        GSMComponents = np.asarray([hp.read_map(s.GSMlocation + "gsm" + str(i+1) + ".fits" + str(s.GSMNSIDE),verbose=False,) for i in range(3)])
        components = np.loadtxt(s.GSMlocation + "components.dat")
        temp = 10**(interp1d(np.log10(components[:,0]), np.log10(components[:,1]), kind='cubic')(np.log10(s.freq))) #cubic spline in log(T)
        weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='linear')(s.freq) for i in range(3)]) #linear interpolation for weights
        GSM = temp*np.dot(weights,GSMComponents)
        interpoltedGSMRotated = hp.get_interp_val(GSM,-g.galCoords.b.radian+np.pi/2, np.asarray(g.galCoords.l.radian))
        hp.mollview(np.log10(interpoltedGSMRotated), title="GSM in Rotated Equatorial Coordinates")


#TODO: investigate temperature/flux normalization        
        
        
#TODO: FIX WHAT HAPPENS WHEN GSM AND MAP RESOLUTION AREN'T THE SAME RESOLUTION  
        
        for t in range(len(s.LSTs)):
            pixelAlts, pixelAzs = mmh.convertEquatorialToHorizontal(s,g.pixelRAs,g.pixelDecs,s.LSTs[t])
            rHatVectors = mmh.convertAltAzToCartesian(pixelAlts,pixelAzs)
            primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[0]), np.pi/2-pixelAlts, pixelAzs)
            for b in range(len(s.baselines)):
                exponent = np.exp(-2j*np.pi*np.dot(rHatVectors,s.baselines[b]))
                visibilities[t,b] += np.sum(GSM * primaryBeam * exponent) * 4*np.pi / len(GSM)
                
#    if s.simulateVisibilitiesWithPointSources:
        #Load PS catalog
        

    return visibilities