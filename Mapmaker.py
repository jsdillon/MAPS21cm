import numpy as np
import healpy as hp
import os
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
import Source.MapmakerHelper as mm


freq = 150
#TODO: the frequency should be passed as a command line parameter

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
s = Specifications(scriptDirectory, "/configuration.txt",freq)
PBs = PrimaryBeams(s)

#other info we need:
	#list of LSTs
	#perVisibilityNoise in jansky
	#field specifications

#print len(s.LSTs)
#mm.CutOutUnusedLSTs(s)
#print len(s.LSTs)

if s.simulateVisibilitiesWithGSM or s.simulateVisibilitiesWithPointSources:
	visibilities = VisibilitySimulator(s,PBs)

from scipy.interpolate import interp1d
import ephem

#Compute GSM from 3 principal components appropriately weighted
GSMComponents = np.asarray([hp.read_map(s.GSMlocation + "gsm" + str(i+1) + ".fits" + str(s.GSMnSide),verbose=False,) for i in range(3)])
components = np.loadtxt(s.GSMlocation + "components.dat")
temp = 10**(interp1d(components[:,0], np.log10(components[:,1]), kind='cubic')(s.freq)) #cubic spline in log(T)
weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='linear')(s.freq) for i in range(3)]) #linear interpolation for weights
GSM = temp*np.dot(weights,GSMComponents)

#Convert GSM to equatorial coordinates
GSMequatorial = GSM

#Let's start by finding the GSM in RA and Dec. Then we'll worry about rotating it properly.

#equaCoords = np.transpose(np.asarray(hp.pix2ang(s.GSMnSide,np.arange(12*s.GSMnSide**2))))

healpixCoords = np.transpose(np.asarray(hp.pix2ang(s.GSMnSide,np.arange(12*s.GSMnSide**2))))
healpixCoords[:,0] = -healpixCoords[:,0] + np.pi/2

from astropy.coordinates import SkyCoord
from astropy import units as u
thisCoord = SkyCoord(frame="galactic", l=healpixCoords[:,1]*u.rad, b=healpixCoords[:,0]*u.rad)
hp.mollview(thisCoord.transform_to("icrs").ra.radian, title="ras")

#hp.mollview(hpCoords[:,1], title="Mollview image RING")
plt.show()


#Find the equatorial coordinates of all points on the equapoint sphere





#hp.mollview(allDecs, title="Mollview image RING")
#plt.show()

#print len(equaCoords)

# Converts a pointing on the healpix sphere (which as been rotated so that the facet is centered on the equator) to an equatorial pointing
# equaPoint convertToEquaPoint(pointing& hpPoint){
# 	double RA = fmod(hpPoint.phi + facetRAinRad - pi, 2*pi);
# 	double dec = pi/2 - hpPoint.theta + facetDecInRad;
# 	if ((dec < -pi/2) || dec > pi/2) RA = fmod(RA + pi, 2*pi);
# 	if (dec < -pi/2) dec = -pi/2 + (-pi/2 - dec);
# 	if (dec > pi/2) dec = pi/2 - (dec - pi/2);
# 	return equaPoint(RA,dec);
# }


#[pixel number, theta/phi]

#print equaCoords[0:5,:]



#hp.mollview(GSM, title="Mollview image RING")


#functions I need:
#	calculate the RA and Dec of each pixel on the sphere where the facet center is placed at the equator and the 
#	Convert those coordinates into galactic coordinates
#	Interpolate on the galactic GSM





#print np.shape(equaCoords[0,:])

