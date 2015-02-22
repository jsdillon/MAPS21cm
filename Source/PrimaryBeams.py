# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp

#When constructed, this class loads in information about the primary beams, using a specifications object
class PrimaryBeams:

    def __init__(self,s):
        print "Now loading in information about the primary beams..."
        self.allBeams = {}
        for antPol in s.antPolList:
            for skyPol in s.skyPolList:
                for antIndex in range((s.nAntennas-1)*(not s.antennasHaveIdenticalBeams) + 1): 
                    #either range(1) or range(nAntennas)
                    for pointIndex in range(s.nPointings):
                        #find two closest values
                        allFrequencies = np.asarray(map(float,s.beamFreqList))
                        freq1Index = np.abs(allFrequencies - s.freq).argmin()
                        freq1 = allFrequencies[freq1Index]
                        allFrequencies[freq1Index] = 1e10
                        freq2Index = np.abs(allFrequencies - s.freq).argmin()
                        freq2 = allFrequencies[freq2Index]
                        
                        filename = s.beamFileFormat.replace('[antIndex]',str(antIndex)).replace('[antPol]',antPol).replace('[skyPol]',skyPol).replace('[pointIndex]',str(pointIndex))
                        key = str(antIndex) + ";" + str(antPol) + ";" + str(skyPol) + ";" + str(pointIndex)
                        beam1 = np.load(filename.replace('[freq]',"{:1.6f}".format(float(s.beamFreqList[freq1Index]))));
                        beam2 = np.load(filename.replace('[freq]',"{:1.6f}".format(float(s.beamFreqList[freq2Index]))));
                        #linear interpolation
                        self.allBeams[key] = beam1 * (1 - (s.freq - freq1)/(freq2 - freq1)) + beam2 * ((s.freq - freq1)/(freq2 - freq1))
        
        #Beam products useful for Stokes I
        for antPol in s.antPolList:
            for skyPol in s.skyPolList:
                for pointIndex in range(s.nPointings):
                    keyIn = str(0) + ";" + str(antPol) + ";" + str(skyPol) + ";" + str(pointIndex)
                    keyOut = str(antPol) + str(antPol) + str(skyPol) + str(skyPol) + str(pointIndex)
                    self.allBeams[keyOut] = self.allBeams[keyIn] * self.allBeams[keyIn].conj()

    #This function returns beam XX, YY, etc. as a numpy array representing a healpix map
    def beamSq(self,antPol,skyPol,point = 0):	
        try:
            return self.allBeams[str(antPol) + str(antPol) + str(skyPol) + str(skyPol) + str(point)]
        except:
            print "Error: cannot find beam " + str(antPol) + str(antPol) + str(skyPol) + str(skyPol) + str(point)
    
    #This function returns beam XX, YY, YX, etc. as a numpy array representing a healpix map    	
    def beamProduct(self,ant1,ant2,antPol1,antPol2,skyPol1,skyPol2,point1 = 0, point2 = 0):
        try:
            key1 = str(ant1) + ";" + str(antPol1) + ";" + str(skyPol1) + ";" + str(point1)
            key2 = str(ant1) + ";" + str(antPol1) + ";" + str(skyPol1) + ";" + str(point2)
            return self.allBeams[key1] * self.allBeam[key2].conj()
        except:
            print "Error: cannot find beam product " + str(ant1) + str(antPol1) + str(skyPol1) + str(point1) + " times " + str(ant1) + str(antPol1) + str(skyPol1) + str(point2) 

