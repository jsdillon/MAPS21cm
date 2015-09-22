# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp

print "Now loading component_maps_408locked.dat"
GSMComponents = np.loadtxt("component_maps_408locked.dat")
for NSIDE in 2**np.arange(3,10):
    print "Working on " + str(NSIDE)
    for comp in range(3):
    	GSMComponentDegraded = hp.pixelfunc.ud_grade(GSMComponents[:,comp],NSIDE)
    	np.save("component_maps_408locked_NSIDE-" + str(NSIDE) + "_Comp-" + str(comp),GSMComponentDegraded)