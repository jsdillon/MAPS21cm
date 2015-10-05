# SUPPORTING FUNCTIONALITY FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import Specifications
import Geometry

def LoadVisibilities(s,times):
    """For now, this just returns a matrix of zeros of the correct size."""
    print "Visibility loading functions are not done. All visibilities set to 0."    
    return np.zeros((len(times.LSTs), s.nBaselines), dtype=np.complex128)
