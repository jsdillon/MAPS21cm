# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import math
import Geometry

class PointSourceCatalog:
    def __init__(self,s,PBs):
        #Compute GSM from 3 principal components appropriately weighted        
        self.freq = s.freq
        self.catalog = np.loadtxt(s.pointSourceCatalogFilename)

        # Determines if the beam-weighted flux of each point source is above the limit and deletes it from the catalog if it isn't
        middleLSTindex = int(math.floor(len(s.LSTs)/2.0))
        psAlts, psAzs = Geometry.convertEquatorialToHorizontal(s, self.catalog[:,0] * 2*np.pi/360, self.catalog[:,1] * 2*np.pi/360, s.LSTs[middleLSTindex])
        primaryBeamWeights = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[middleLSTindex]), np.pi/2-psAlts, psAzs)
        beamWeightedFluxes = primaryBeamWeights * self.catalog[:,2] * (s.freq/s.pointSourceReferenceFreq)**(-self.catalog[:,3])
        self.catalog = self.catalog[beamWeightedFluxes > s.pointSourceBeamWeightedFluxLimit, :]
        print str(len(self.catalog)) + " point sources identified for specific modeling."
        
        #Convert into a more useful format        
        self.RAs = self.catalog[:,0] * 2*np.pi/360
        self.decs = self.catalog[:,1] * 2*np.pi/360
        self.fluxes = self.catalog[:,2] 
        self.spectralIndices = self.catalog[:,3]    
        self.scaledFluxes = self.fluxes * (s.freq / s.pointSourceReferenceFreq)**(-self.spectralIndices)