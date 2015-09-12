# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This file sets up a description of the observation that would be provided earlier in the pipeline for real data

import os
import numpy as np
import ConfigParser
import pickle

#HARD CODED SETTINGS
integrationTime = 10 #seconds
singleBaselineVisibilityNoise = 1.0 #Jy

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
LSTs = np.asarray(np.arange(0,24.0,integrationTime/60.0/60.0))
np.savetxt(scriptDirectory + "/LSTs.dat",LSTs)

config = ConfigParser.ConfigParser()
config.read(scriptDirectory + '/../configuration.txt')
useOnlyUniqueBaselines = config.getboolean('Array Settings','useOnlyUniqueBaselines')
antennasHaveIdenticalBeams = config.getboolean('Array Settings','antennasHaveIdenticalBeams')
mainDirectory = scriptDirectory + '/..'
antennaPositions = np.loadtxt(config.get('Array Settings','antennaPositionsFile').replace('[MainDirectory]',mainDirectory))
allBaselines = np.loadtxt(config.get('Array Settings','allBaselinesListFile').replace('[MainDirectory]',mainDirectory))
uniqueBaselines = np.loadtxt(config.get('Array Settings','uniqueBaselinesListFile').replace('[MainDirectory]',mainDirectory))
baselineRedundancies = np.loadtxt(config.get('Array Settings','baselineRedundancyFile').replace('[MainDirectory]',mainDirectory))


#all beams taken to be pointed to position 0 for now
if antennasHaveIdenticalBeams: #here we assume that identical beams also means identical pointings at any given time
	pointings = np.zeros(np.shape(LSTs))
else:
	pointings = np.zeros((len(LSTs),len(antennaPositions)))
np.savetxt(scriptDirectory + "/pointings.dat",pointings,fmt='%i')

pointingCenters = {0: [np.pi/2, 0]}
pickle.dump(pointingCenters, open('pointing_centers.p', 'w'))

#all beams are taken to be pointed at the zenith (alt = pi/2, az = 0)

# If we're binning together redundant baselines, then we want to know the noise on that baseline during the integration
if useOnlyUniqueBaselines:
	baselineNoise = np.ones((len(LSTs),len(uniqueBaselines)))*singleBaselineVisibilityNoise / (baselineRedundancies**.5)
	np.save(scriptDirectory + "/noisePerUniqueBaseline",baselineNoise)	
# If we're not binning together redundant baselines, then we want to know the per-antenna noise and then take the product of the two antennas in the baseline
else:
	antennaNoise = np.ones((len(LSTs),len(antennaPositions))) * (singleBaselineVisibilityNoise**.5)
	np.save(scriptDirectory + "/noisePerAntenna", antennaNoise)


#TODO: 
#	HERA will give me all the baselines and let me sum them
#	HERA will give me weights (integers) as a funciton of time, frequency, and baseline
#	I'll have to calculate antenna noise as a function of (LS) time
# 	Assumes no other variation in noise temperature
#	Weights cube will be a uv FITS like everything else