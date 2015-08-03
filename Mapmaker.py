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
    #visibilties = mmh.LoadVisibilities(s)
visibilities *= s.convertJyToKFactor

Geometry.rephaseVisibilitiesToSnapshotCenter(s,visibilities,times)
MapMats.inverseCovarianceWeightVisibilities(s,visibilities)

for snapshot in times.snapshots:    
    NinvTimesy = MapMats.calculateNinvTimesy(visibilities, snapshot)
    Ninv = MapMats.calculateNInv(s,snapshot)
    
}    
    
#	Healpix_Map<double> coaddedMap = emptyHealpixMap();
#	vector< vector<double> > PSF(nPixels, vector<double>(nPixelsExtended,0.0));
#	vector< vector<double> > pointSourcePSF(nPixels, vector<double>(nPointSources,0.0));
#	for (int n = 0; n < nSnapshots; n++){
#		cout << " " << floor(100.0 * n / nSnapshots) << "% done. \r" << std::flush;
#		int snapshotCentralLSTindex = snapshotLSTindices[n][int(round(snapshotLSTindices[n].size()/2.0-.5))];
#		vector<complex> NinvTimesy = calculateNinvTimesy(allVisibilities, snapshotLSTindices[n]);
#		vector<double> Ninv = calculateNinv(noiseVarianceOnEachVisibiltiy, snapshotLSTindices[n]);
#		vector< vector<complex> > KAtranspose = calculateKAtranspose(LSTs[snapshotCentralLSTindex], extendedPixelEquaPointings, baselines, discretizedPrimaryBeam);
#		addSnapshotMap(coaddedMap, NinvTimesy, KAtranspose, mapOfIndicesInExtendedIndexVector, healpixIndices);
#		addSnapshotPSF(PSF, KAtranspose, Ninv, mapOfIndicesInExtendedIndexVector);
#		if (alsoComputePointSourcePSF){
#			vector< vector<complex> > pointSourceAmatrix = calcualtePointSourceAmatrix(LSTs[snapshotCentralLSTindex], baselines, discretizedPrimaryBeam, allPointSources);
#			addSnapshotPointSourcePSF(pointSourcePSF, KAtranspose, Ninv, pointSourceAmatrix, mapOfIndicesInExtendedIndexVector);
#		}
#	}



#PSUEDO CODE:
#-loop over snapshots:
#    -calculate Ninv DONE
#    -Calculate Ninv*y
#    -Calculate KAtranspose
#    -add snapshot map
#    -add snapshot PSF
#Normalize the final output
#Save data products    


# Packages needed: Mapmaking Matrices, 

#plt.figure()
#plt.plot(np.real(visibilities))



#hp.mollview(pixelAzs, title="Pixel Azimuths at LST=0")
#hp.mollview(PBs.beamSquared("X","x",s.pointings[0]), title="PrimaryBeam")

plt.show()    

    
