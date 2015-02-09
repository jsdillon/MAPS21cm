# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import ConfigParser
import numpy as np

#This class takes the location of the configuration file and loads relevant specs as attributes to the object
class Specifications:
    def __init__(self,directory,configFilename,freq):
        
        self.freq = freq
        config = ConfigParser.ConfigParser()
        self.mainDirectory = directory
        config.read(directory + configFilename)

        try: 
			self.useOnlyUniqueBaselines = config.getboolean('Array Settings','useOnlyUniqueBaselines')
			self.antennaPositions = np.loadtxt(config.get('Array Settings','antennaPositionsFile').replace('[MainDirectory]',self.mainDirectory))
			self.allBaselines = np.loadtxt(config.get('Array Settings','allBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
			self.uniqueBaselines = np.loadtxt(config.get('Array Settings','uniqueBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
			self.baselineRedundancies = np.loadtxt(config.get('Array Settings','baselineRedundancyFile').replace('[MainDirectory]',self.mainDirectory))

			self.antennasHaveIdenticalBeams = config.getboolean('Array Settings','antennasHaveIdenticalBeams')
			self.antPolList = config.get('Array Settings','antPolList').split()
			self.skyPolList = config.get('Array Settings','skyPolList').split()
			self.nPointings = config.getint('Array Settings','nPointings')
			self.beamFreqList = config.get('Array Settings','beamFreqList').split()
			self.beamFileFormat = config.get('Array Settings','beamFileFormat').replace('[MainDirectory]',self.mainDirectory)
			self.beamNSIDE = config.getint('Array Settings','beamNSIDE')
			
			self.simulateVisibilities = config.getboolean('Input Data Settings','simulateVisibilities')

        except:
            print "\nWARNING: Something is wrong with configuration.txt\n"

        #Other calculations based on inputs
        self.nAntennas = len(self.antennaPositions)


        
