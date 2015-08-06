# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import math
import Geometry
import cPickle as pickle
import os

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

def calculateKAtranspose(s,snapshot,coords,PBs):
    """This function computes K_PSF * A^t, which maps baselines at the given snapshot central index to the extended facet.\n
    It is worth nothing that the Fourier convention is based on A having e^ib.k, so A^t has e^-ib.k = e^i|k|b.theta_hat, as we see here."""
    realSpaceDiagonalPart = np.ones(coords.nExtendedPixels) * 4*np.pi / 12.0 / coords.NSIDE**2 * s.convertJyToKFactor
    extendedAlts, extendedAzs = Geometry.convertEquatorialToHorizontal(s, coords.pixelRAs[coords.extendedIndices], coords.pixelDecs[coords.extendedIndices],snapshot.centralLST)
    extendedCartVecs = Geometry.convertAltAzToCartesian(extendedAlts, extendedAzs)
    realSpaceDiagonalPart *= hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[snapshot.centralLSTIndex]), np.pi/2-extendedAlts, extendedAzs)    
    KAtranspose = np.dot(np.diag(realSpaceDiagonalPart), np.exp(1j * s.k * extendedCartVecs.dot(np.transpose(s.baselines))))
    return KAtranspose
    
def calculatePSAmatrix(s,snapshot,ps,PBs):
    """This function computes A mappings to the locations of the point source at the snapshot central index to the baselines."""
    psAlts, psAzs = Geometry.convertEquatorialToHorizontal(s, ps.RAs, ps.decs, snapshot.centralLST)
    psCartVecs = Geometry.convertAltAzToCartesian(psAlts, psAzs)
    realSpaceDiagonalPart = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[snapshot.centralLSTIndex]), np.pi/2-psAlts, psAzs)
    pointSourceAmatrix = np.dot(np.exp(-1j * s.k * s.baselines.dot(np.transpose(psCartVecs))), np.diag(realSpaceDiagonalPart))
    return pointSourceAmatrix
    
def saveAllResults(s,coords,times,ps,Dmatrix,PSF,coaddedMap,pointSourcePSF,mapNoiseCovariance):
    """This function saves all the input classes, vectors, and matrices to s.resultsFolder."""
    os.system("rm -rf " + s.resultsFolder)
    os.system("mkdir " + s.resultsFolder)
    pickle.dump(s, open(s.resultsFolder + "specifications.p","wb"))
    pickle.dump(coords, open(s.resultsFolder + "coordinates.p","wb"))
    pickle.dump(times, open(s.resultsFolder + "times.p","wb"))
    pickle.dump(ps, open(s.resultsFolder + "pointSourceCatalog.p","wb"))
    np.save(s.resultsFolder + "Dmatrix",Dmatrix)
    np.save(s.resultsFolder + "PSF",PSF)
    np.save(s.resultsFolder + "coaddedMap",coaddedMap)
    np.save(s.resultsFolder + "pointSourcePSF",pointSourcePSF)
    np.save(s.resultsFolder + "mapNoiseCovariance",mapNoiseCovariance)

