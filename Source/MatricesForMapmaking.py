# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import math

def inverseCovarianceWeightVisibilities(s,visibilities):
    """This function weights the visibilities (which have not been combined into snapshots yet) by the inverse of the noise variance."""
    if s.useOnlyUniqueBaselines:
        visibilities /= (s.noisePerUniqueBaseline**2)
    else:
        print "WARNING: this has not been tested."        
        for b in range(len(s.baselines)):
            visibilities[:,b] /= (s.noisePerAntenna[:,s.allBaselinePairs[b,0]] * s.noisePerAntenna[:,s.allBaselinePairs[b,1]])

def calculateNinvTimesy(visibilities, snapshot):
    """This function computes the inverse varianced weighted sum of rephased visibilities in a snapshot. The averaging (rather than summing) happens when the normalizaiton is applied."""
    return np.sum(visibilities[snapshot.LSTindices,:], axis=0)
    
def calculateNInv(s, snapshot):
    """This function computes the diagonal of the inverse noise covariance (off-diagonal terms are 0) for each baseline during this snapshot."""
    if s.useOnlyUniqueBaselines:
        return np.sum(s.noisePerUniqueBaseline[snapshot.LSTindices,:]**(-2), axis=0)
    else:
        print "WARNING: this has not been tested."
        return np.asarray([np.sum((s.noisePerAntenna[:,s.allBaselinePairs[b,0]] * s.noisePerAntenna[:,s.allBaselinePairs[b,1]])**(-1), axis=0) for b in range(len(s.baselines))])
