# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import ConfigParser
import numpy as np
import math
import ephem
import pickle
import scipy.constants as const

#This class takes the location of the configuration file and loads relevant specs as attributes to the object
class Specifications:
    def __init__(self,directory,configFilename,freq):
        print "Now loading in specifications..."
        self.freq = freq
        self.k = 2*np.pi * freq*1e6 / const.c
        config = ConfigParser.ConfigParser()
        self.mainDirectory = directory
        config.read(directory + configFilename)

        #ANTENNA AND INSTRUMENT SETTINGS
        self.useOnlyUniqueBaselines = config.getboolean('Array Settings','useOnlyUniqueBaselines')
        self.antennaPositions = np.loadtxt(config.get('Array Settings','antennaPositionsFile').replace('[MainDirectory]',self.mainDirectory))
        if self.useOnlyUniqueBaselines:
            self.baselines = np.loadtxt(config.get('Array Settings','uniqueBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
            self.baselineRedundancies = np.loadtxt(config.get('Array Settings','baselineRedundancyFile').replace('[MainDirectory]',self.mainDirectory))
        else:
            self.baselines = np.loadtxt(config.get('Array Settings','allBaselinesListFile').replace('[MainDirectory]',self.mainDirectory))
            self.allBaselinePairs = np.loadtxt(config.get('Array Settings','allBaselinePairsListFile').replace('[MainDirectory]',self.mainDirectory))
            self.baselineRedundancies = np.ones(self.baselines.shape[0]) #since no baselines are redundant
        self.arrayLat = config.getfloat('Array Settings','arrayLat')
        self.arrayLatInRad = self.arrayLat * math.pi/180.0;
        self.arrayLong = config.getfloat('Array Settings','arrayLong')

        
        #BEAM SETTINGS
        self.antennasHaveIdenticalBeams = config.getboolean('Array Settings','antennasHaveIdenticalBeams')
        if (self.useOnlyUniqueBaselines and not self.antennasHaveIdenticalBeams):
            print "\nWARNING: Cannot use group baselines identically if the beams are not identical.\n"
        self.antPolList = config.get('Array Settings','antPolList').split()
        self.skyPolList = config.get('Array Settings','skyPolList').split()        
        self.beamFreqList = config.get('Array Settings','beamFreqList').split()
        self.beamFileFormat = config.get('Array Settings','beamFileFormat').replace('[MainDirectory]',self.mainDirectory)
        self.beamNSIDE = config.getint('Array Settings','beamNSIDE')    
        
        #OBSERVATION SETTINGS
        self.LSTsFilename = config.get('Input Data Settings','LSTsFilename').replace('[MainDirectory]',self.mainDirectory)
        self.pointings = np.loadtxt(config.get('Input Data Settings','PointingListFilename').replace('[MainDirectory]',self.mainDirectory)).astype(int)
        self.pointingCenters = pickle.load(open(config.get('Input Data Settings','PointingCenterDictionaryFilename').replace('[MainDirectory]',self.mainDirectory),'r'))
        if self.useOnlyUniqueBaselines:
            self.noisePerUniqueBaseline = np.load(config.get('Input Data Settings','noisePerUniqueBaselineFilename').replace('[MainDirectory]',self.mainDirectory))
        else:
            self.noisePerAntenna = np.load(config.get('Input Data Settings','noisePerAntennaFilename').replace('[MainDirectory]',self.mainDirectory))

        #POINT SOURCE CATALOG SETTINGS
        self.pointSourceCatalogFilename = config.get('Input Data Settings','pointSourceCatalogFilename').replace('[MainDirectory]',self.mainDirectory)
        self.pointSourceReferenceFreq = config.getfloat('Input Data Settings','pointSourceReferenceFreq')
        self.pointSourceBeamWeightedFluxLimitAtReferenceFreq = config.getfloat('Input Data Settings','pointSourceBeamWeightedFluxLimitAtReferenceFreq')
        
        
        #VISIBILITY SIMULATION SETTINGS
        self.simulateVisibilitiesWithGSM = config.getboolean('Input Data Settings','simulateVisibilitiesWithGSM')
        self.simulateVisibilitiesWithPointSources = config.getboolean('Input Data Settings','simulateVisibilitiesWithPointSources')
        self.GSMlocation = config.get('Input Data Settings','GSMlocation').replace('[MainDirectory]',self.mainDirectory)
        self.GSMNSIDE = config.getint('Input Data Settings','GSMNSIDE')
        
        #FACET SETTINGS
        self.facetRA = config.getfloat('Mapmaking Specifications','facetRA')
        self.facetRAinRad = self.facetRA * math.pi/180.0;
        self.facetDec = config.getfloat('Mapmaking Specifications','facetDec')
        self.facetDecinRad = self.facetDec * math.pi/180.0;
        self.facetSize = config.getfloat('Mapmaking Specifications','facetSize')
        self.MaximumAllowedAngleFromFacetCenterToPointingCenter = config.getfloat('Mapmaking Specifications','MaximumAllowedAngleFromFacetCenterToPointingCenter')
        
        #MAPMAKING AND PSF SETTINGS        
        self.mapNSIDE = config.getint('Mapmaking Specifications', 'mapNSIDE')
        self.PSFextensionBeyondFacetFactor = config.getfloat('Mapmaking Specifications', 'PSFextensionBeyondFacetFactor')    
        self.integrationsPerSnapshot = config.getint('Mapmaking Specifications', 'integrationsPerSnapshot')
        self.PSFforPointSources = config.getboolean('Mapmaking Specifications','PSFforPointSources')
        self.mapPixels = 12 * self.mapNSIDE**2
        
        #OUTPUT SETTINGS
        self.resultsFolder = config.get('Mapmaking Specifications','resultsFolder').replace('[MainDirectory]',self.mainDirectory) + "{:.3f}".format(self.freq) + "/"        

        #Other calculations based on inputs
        self.nAntennas = len(self.antennaPositions)
        self.nPointings = len(self.pointingCenters)
        self.convertJyToKFactor = (const.c)**2 / (2 * const.k * (self.freq*1e6)**2 * 1e26) #multiply by this number to convert Jy to K