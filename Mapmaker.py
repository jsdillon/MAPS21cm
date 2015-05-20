import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
import Source.MapmakerHelper as mmh
import scipy.constants as const
plt.close("all")

freq = 150
#TODO: the frequency should be passed as a command line parameter

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
s = Specifications(scriptDirectory, "/configuration.txt",freq)
geo = mmh.Geometry(s)
PBs = PrimaryBeams(s)

mmh.CutOutUnusedLSTs(s)
mmh.DecideWhichSourcesToInclude(s,PBs)




#other info we need:	#list of LSTs
	#perVisibilityNoise in jansky
	#field specifications

if s.simulateVisibilitiesWithGSM or s.simulateVisibilitiesWithPointSources:
    visibilities = VisibilitySimulator(s,PBs)
else:
    visibilties = mmh.LoadVisibilities(s)
visibilities *= s.convertJyToKFactor

plt.figure()
plt.plot(np.real(visibilities))

#TODO: CONVERT VISIBILITIES INTO TEMPERATURE UNITS


#hp.mollview(pixelAzs, title="Pixel Azimuths at LST=0")

#Visibility simulation

#hp.mollview(PBs.beamSquared("X","x",s.pointings[0]), title="PrimaryBeam")

plt.show()    

    
