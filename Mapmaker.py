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
GSMComponents = np.asarray([hp.read_map(s.GSMlocation + "gsm" + str(i+1) + ".fits" + str(s.GSMNSIDE),verbose=False,) for i in range(3)])
components = np.loadtxt(s.GSMlocation + "components.dat")
temp = 10**(interp1d(components[:,0], np.log10(components[:,1]), kind='cubic')(s.freq)) #cubic spline in log(T)
weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='linear')(s.freq) for i in range(3)]) #linear interpolation for weights
GSM = temp*np.dot(weights,GSMComponents)



from astropy.coordinates import SkyCoord
from astropy import units as u

healpixCoords = np.transpose(np.asarray(hp.pix2ang(s.mapNSIDE,np.arange(s.mapPixels))))

#What I really want is to take a set of healpix vectors where (x=1,y=0,z=0) is the facet center. 
	# This is really just
	# I want to apply a rotation matrix 
	# 

#No rotation to facet center
pixelRAs = healpixCoords[:,1]
pixelDecs = -healpixCoords[:,0] + np.pi/2
pixelVectors = [ np.array([np.cos(pixelDecs[n]) * np.cos(pixelRAs[n]), np.cos(pixelDecs[n]) * np.sin(pixelRAs[n]), np.sin(pixelDecs[n])]) for n in range(12*s.GSMNSIDE**2)]

rotateToPrimeMeridian = np.array([[np.cos(s.facetRAinRad), -np.sin(s.facetRAinRad), 0], [np.sin(s.facetRAinRad), np.cos(s.facetRAinRad), 0],[0, 0, 1]])
rotateToEquator = np.array([[np.cos(-s.facetDecinRad), 0, np.sin(-s.facetDecinRad)], [0, 1, 0], [-np.sin(-s.facetDecinRad), 0, np.cos(-s.facetDecinRad)]])
#rotatedVectors = np.asarray([np.dot(np.dot(rotateToEquator,rotateToPrimeMeridian),pixelVectors[n]) for n in range(s.mapPixels)])
rotatedVectors = [np.dot(np.dot(rotateToEquator,rotateToPrimeMeridian),pixelVectors[n]) for n in range(s.mapPixels)]

	

newPixelDecs = np.asarray([np.arcsin(rotatedVectors[n][2]) for n in range(s.mapPixels)])
newPixelRAs = np.asarray([np.arctan2(rotatedVectors[n][1],rotatedVectors[n][0]) for n in range(s.mapPixels)])
newPixelRAs[newPixelRAs < 0] = newPixelRAs[newPixelRAs < 0] + 2*np.pi

hp.mollview(pixelRAs,title="pixelRAs")
hp.mollview(newPixelRAs, title="newPixelRAs")
#hp.mollview(newPixelDecs, title="newPixelDecs")
#Next rotate so that facet center is on the equator

#hp.mollview(pixelRAs, title="Pixel RAs")
#hp.mollview(rotatedVectors[:,0], title="Pixel xs")
#hp.mollview(rotatedVectors[:,1], title="Pixel ys")
#hp.mollview(rotatedVectors[:,2], title="Pixel zs")
#hp.mollview(pixelRAs, title="Pixel RAs")
#hp.mollview(pixelYs, title="Pixel Decs")

#print 



galCoords = SkyCoord(frame="icrs", ra=pixelRAs*u.rad, dec=pixelDecs*u.rad).transform_to("galactic")
interpoltedGSM = hp.get_interp_val(GSM,-galCoords.b.radian+np.pi/2, np.asarray(galCoords.l.radian))
#hp.mollview(np.log10(interpoltedGSM), title="GSM in Equatorial Coordinates")



plt.show()

