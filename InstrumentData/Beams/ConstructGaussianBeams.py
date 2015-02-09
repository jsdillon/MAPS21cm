# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This code creates a gaussian antenna beam with hard-coded FWHM at 150 MHz (of the beam squared) pointed at the zenith

import numpy as np
import healpy as hp
import math

NSIDE = 64
freqList = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
FWHMat150 = 10
antennaIndex = 0
feedPolarizationList = ["X","Y"]
skyPolarizationlist = ["x","y"]
pointIndex = 0

for freq in freqList:
	print "Now working on " + str(freq) + " MHz..."
	for pol in feedPolarizationList:
		for skyPol in skyPolarizationlist:
			FWHM = FWHMat150 * (150.0/freq)
			sigma = (FWHM/2.35482004503) * (math.pi/180)
			antennaBeamValues = np.zeros(12*NSIDE*NSIDE)

			for n in range(12*NSIDE*NSIDE):
				angleHere = hp.pix2ang(NSIDE,n)
				if angleHere[0] < math.pi/2:
					angularDistance = hp.rotator.angdist(angleHere,[0,0])[0]
					antennaBeamValues[n] = math.exp(-angularDistance**2 / (2*2 * sigma**2)) #This is so the antenna beam squared has the right FWMH

			outfilename = "beam_" + str(antennaIndex) + "_" + pol + "_" + skyPol + "_" +  str(pointIndex) + "_{:1.6f}".format(freq)
			np.save(outfilename,antennaBeamValues)

print "Done!\n"