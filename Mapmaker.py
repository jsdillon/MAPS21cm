import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
import scipy.constants as const
plt.close("all")

freq = 150
#TODO: the frequency should be passed as a command line parameter

#Load in everything we need, figure out which LSTs to work with
scriptDirectory = os.path.dirname(os.path.abspath(__file__))
s = Specifications(scriptDirectory, "/configuration.txt",freq)
times = Geometry.Times(s)
times.CutOutUnusedLSTsAndGroupIntoSnapshots(s)
coords = Geometry.Coordinates(s)
PBs = PrimaryBeams(s)
ps = PointSourceCatalog(s,PBs,times)


#Simulate or load visibilities
if s.simulateVisibilitiesWithGSM or s.simulateVisibilitiesWithPointSources:
    visibilities = VisibilitySimulator(s,PBs,ps,times)
else:
    print "Visibility loading functions are not done."    
    #visibilties = LoadVisibilities(s)
visibilities *= s.convertJyToKFactor

#Prepare visibilities
Geometry.rephaseVisibilitiesToSnapshotCenter(s,visibilities,times)
MapMats.inverseCovarianceWeightVisibilities(s,visibilities)

#Perform mapmaking and calculate PSFs
coaddedMap = np.zeros(coords.nFacetPixels)
PSF = np.zeros((coords.nFacetPixels,coords.nExtendedPixels))
pointSourcePSF = np.zeros((coords.nFacetPixels, ps.nSources))
for snapshot in times.snapshots:    
    NinvTimesy = MapMats.calculateNinvTimesy(visibilities, snapshot)
    Ninv = MapMats.calculateNInv(s,snapshot)
    KAtranspose = MapMats.calculateKAtranspose(s,snapshot,coords,PBs)    
    coaddedMap += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], NinvTimesy))
    PSF += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], np.dot(np.diag(Ninv), KAtranspose.conj().transpose())))    
    if s.PSFforPointSources and ps.nSources > 0:
        pointSourceAmatrix = MapMats.calculatePSAmatrix(s,snapshot,ps,PBs)
        pointSourcePSF += 2 * np.real(np.dot(KAtranspose[coords.mapIndexLocationsInExtendedIndexList,:], np.dot(np.diag(Ninv), pointSourceAmatrix)))
        
#Renormalize maps and PSFs
Dmatrix = np.diag(np.diag(PSF[:,coords.mapIndexLocationsInExtendedIndexList])**(-1))
PSF = np.dot(Dmatrix,PSF)
coaddedMap = np.dot(Dmatrix,coaddedMap)
pointSourcePSF = np.dot(Dmatrix,pointSourcePSF)
mapNoiseCovariance = np.dot(PSF[:,coords.mapIndexLocationsInExtendedIndexList],np.transpose(Dmatrix))
 
    




#PSUEDO CODE:
#Save data products    



plt.show()    

    
