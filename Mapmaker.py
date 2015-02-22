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
GSMComponents = np.asarray([hp.read_map(s.GSMlocation + "gsm" + str(i+1) + ".fits" + str(s.GSMnSide)) for i in range(3)])
components = np.loadtxt(s.GSMlocation + "components.dat")
temp = 10**(interp1d(components[:,0], np.log10(components[:,1]), kind='cubic')(s.freq)) #cubic spline in log(T)
weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='linear')(s.freq) for i in range(3)]) #linear interpolation for weights
GSM = temp*np.dot(weights,GSMComponents)

#Convert GSM to equatorial coordinates
GSMequatorial = GSM
equaCoords = hp.pix2ang(s.GSMnSide,np.arange(12*s.GSMnSide**2))
plt.close(1)
#print np.shape(equaCoords[0,:])

