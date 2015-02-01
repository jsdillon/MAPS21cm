# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp

#When constructed, this class loads in information about the primary beams, using a specifications object
class PrimaryBeams:
    def __init__(self,s,freq):
        allBeams
        #TODO: compute and save all single antenna beams at this frequency and all XX, YY, etc.

	#This function returns beam XX, YY, etc. as a numpy array representing a healpix map
    def beam(self,pol,point = 0):	
    	try:
    		return allBeams[str(pol)+";"+str(point)]
    	except:
    		print "Error: cannot find beam " + str(pol)+";"+str(point)
    
    #This function returns beam XX, YY, etc. as a numpy array representing a healpix map    	
    def beamProduct(self,ant1,ant2,pol1,pol2,point1 = 0, point2 = 0):
		try:
    		#TODO: calculate and return beam product
    	except:
    		print "Error: cannot find beam product " + str(ant1) + str(pol1) + str(point1) + " times " + str(ant2) + str(pol2) + str(point2)
	
