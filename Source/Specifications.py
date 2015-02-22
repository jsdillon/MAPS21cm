# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import ConfigParser
import numpy as np
import math
import ephem

#This class takes the location of the configuration file and loads relevant specs as attributes to the object
class Specifications:
    def __init__(self,directory,configFilename,freq):
        print "Now loading in specifications..."
        self.freq = freq
        config = ConfigParser.ConfigParser()
        self.mainDirectory = directory
        config.read(directory + configFilename)

        #ANTENNA AND INSTRUMENT SETTINGS
        self.useOnlyUniqueBaselines = config.getboolean('Array Settings','useOnlyUniqueBaselines')
        self.antennaPositions = np.loadtxt(config.get('Array Settings','antennaPositionsFile').replace('[MainDirectory]',self.mainDirectory))
        self.allBaselines = np.loadtxt(config.get('Array Settings','allBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
        self.uniqueBaselines = np.loadtxt(config.get('Array Settings','uniqueBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
        self.baselineRedundancies = np.loadtxt(config.get('Array Settings','baselineRedundancyFile').replace('[MainDirectory]',self.mainDirectory))
        self.arrayLat = config.getfloat('Array Settings','arrayLat')
        self.arrayLatInRad = self.arrayLat * math.pi/180.0;
        self.arrayLong = config.getfloat('Array Settings','arrayLong')

        
        #BEAM SETTINGS
        self.antennasHaveIdenticalBeams = config.getboolean('Array Settings','antennasHaveIdenticalBeams')
        if (self.useOnlyUniqueBaselines and not self.antennasHaveIdenticalBeams):
            print "\nWARNING: Cannot use group baselines identically if the beams are not identical.\n"
        self.antPolList = config.get('Array Settings','antPolList').split()
        self.skyPolList = config.get('Array Settings','skyPolList').split()
        self.nPointings = config.getint('Array Settings','nPointings')
        self.beamFreqList = config.get('Array Settings','beamFreqList').split()
        self.beamFileFormat = config.get('Array Settings','beamFileFormat').replace('[MainDirectory]',self.mainDirectory)
        self.beamNSIDE = config.getint('Array Settings','beamNSIDE')    
        
        #OBSERVATION SETTINGS
        self.LSTs = np.loadtxt(config.get('Input Data Settings','LSTsFilename').replace('[MainDirectory]',self.mainDirectory))
        self.pointings = np.loadtxt(config.get('Input Data Settings','PointingListFilename').replace('[MainDirectory]',self.mainDirectory))
        self.noiseList = np.load(config.get('Input Data Settings','noiseFilename').replace('[MainDirectory]',self.mainDirectory))

        #VISIBILITY SIMULATION SETTINGS
        self.simulateVisibilitiesWithGSM = config.getboolean('Input Data Settings','simulateVisibilitiesWithGSM')
        self.simulateVisibilitiesWithPointSources = config.getboolean('Input Data Settings','simulateVisibilitiesWithPointSources')
        self.GSMlocation = config.get('Input Data Settings','GSMlocation').replace('[MainDirectory]',self.mainDirectory)
        self.GSMnSide = config.getint('Input Data Settings','GSMnSide')
        self.pointSourceCatalogFilename = config.get('Input Data Settings','pointSourceCatalogFilename').replace('[MainDirectory]',self.mainDirectory)

        #FACET SETTINGS
        self.facetRA = config.getfloat('Mapmaking Specifications','facetRA')
        self.facetRAinRad = self.facetRA * math.pi/180.0;
        self.facetDec = config.getfloat('Mapmaking Specifications','facetDec')
        self.facetDecinRad = self.facetDec * math.pi/180.0;
        self.facetSize = config.getfloat('Mapmaking Specifications','facetSize')
        self.MaximumAllowedAngleFromFacetCenterToZenith = config.getfloat('Mapmaking Specifications','MaximumAllowedAngleFromFacetCenterToZenith')

        #Other calculations based on inputs
        self.nAntennas = len(self.antennaPositions)