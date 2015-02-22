# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import ephem
#from astropy import units as u
#from astropy.coordinates import SkyCoord


def VisibilitySimulator(s,PBs):
    print "testing VisibilitySimulator"


    

    #ephem.Equatorial(equaCoords[0,1],equaCoords[0,2])
    #%c = SkyCoord(equaCoords[0,1], equaCoords[0,2], frame='icrs', unit='rad')
    #print c
    #for n in range(12*s.GSMnSide**2):



    #hp.get_inter_val()


    #hp.mollview(GSM, title="Mollview image RING", norm="log")
    #plt.show()

# //Rebins the GSM into equatorial coordinates, making for an easier comparison between the true sky
# Healpix_Map<double> computeGSMInEquatorialCoordinates(vector<int>& extendedHealpixIndices){
#     int GSMresTemp = GSMres;
#     if (GSMres > 512) GSMres = 512;
#     Healpix_Map<double> galacticGSM = computeGSM();
#     GSMres = GSMresTemp;

#     Healpix_Map<double> equatorialGSM(int(round(log(GSMres)/log(2))), RING);
#     for (int n = 0; n < 12*GSMres*GSMres; n++){
#         pointing thisPoint = equatorialGSM.pix2ang(n);
#         equaPoint thisEquaPoint = convertToEquaPoint(thisPoint);
#         //equaPoint thisEquaPoint(thisPoint.phi, pi/2-thisPoint.theta);
#         pointing thisGalPoint = convertEquaPointToGal(thisEquaPoint);
#         equatorialGSM[n] = galacticGSM.interpolated_value(thisGalPoint);
#     }

#     double sum = 0;
#     for (int n = 0; n < 12*GSMres*GSMres; n++) sum += equatorialGSM[n];
#     cout << "Average temperature after rebinning is " << sum / 12 / GSMres / GSMres << " K." << endl;

#     Healpix_Map<double> galacticGSMlowres = emptyHealpixMap();
#     if (GSMres > NSIDE){
#         galacticGSMlowres.Import_degrade(equatorialGSM);    
#     } else {
#         galacticGSMlowres = equatorialGSM;
#     }
    
#     /*vector<double> wholeGSM(NSIDE*NSIDE*12,0.0);
#     for (int n = 0; n < NSIDE*NSIDE*12; n++) wholeGSM[n] = equatorialGSM[n];
#     exportVector(wholeGSM, NSIDE*NSIDE*12, "wholeGSM.dat"); */


#     vector<double> facetMap(nPixelsExtended, 0.0);
#     for (int n = 0; n < nPixelsExtended; n++) facetMap[n] = galacticGSMlowres[extendedHealpixIndices[n]];
#     exportVector(facetMap, nPixelsExtended, trueSkyFilename);
#     string onlyGSMFilename = dataProductFilePrefix + "trueSkyOnlyGSM.dat";
#     exportVector(facetMap, nPixelsExtended, onlyGSMFilename);
#     return equatorialGSM;
# }





    #Psuedocode:
    #Compute GSM
    #Put GSM into equatorial coordinates
    #Calculate Visibilities
