# TEST CODE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon
#
# This code generates a file of antenna positions and a files with baselines, baseline info, and redundancies

import os
import numpy as np
import cPickle as pickle

#HARD CODED SETTINGS
Separation = 14.6
hexNum = 4
precisionFactor = 1000000

positions = [];
for row in range(hexNum-1,-(hexNum),-1):
	for col in range(2*hexNum-abs(row)-1):
		xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*Separation;
		yPos = row*Separation*np.sqrt(3)/2;
		positions.append([xPos, yPos, 0])

nAntennas = len(positions)
baselines = []
baselinePairs = []
for ant1 in range(nAntennas):
	for ant2 in range(ant1+1,nAntennas):
		baselines.append((int(np.round(precisionFactor*(positions[ant1][0]-positions[ant2][0]))), int(np.round(precisionFactor*(positions[ant1][1]-positions[ant2][1]))), int(np.round(precisionFactor*(positions[ant1][2]-positions[ant2][2])))))
		baselinePairs.append((ant1, ant2))

baselineDict = {}
for b in range(len(baselines)):
	if baselineDict.has_key(baselines[b]):
		baselineDict[baselines[b]].append(baselinePairs[b])
	else:
		baselineDict[baselines[b]] = [baselinePairs[b]]

print "With", len(positions), "antennas there are", len(baselineDict.items()), "unique baselines."

scriptDirectory = os.path.dirname(os.path.abspath(__file__))
np.savetxt(scriptDirectory + "/antenna_positions.dat",np.asarray(positions))
np.savetxt(scriptDirectory + "/all_baselines.dat",np.asarray(baselines)/(1.0*precisionFactor))
np.savetxt(scriptDirectory + "/all_baseline_pairs.dat",np.asarray(baselinePairs),fmt='%i')
np.savetxt(scriptDirectory + "/unique_baselines.dat",np.asarray([uniqueBaseline[0] for uniqueBaseline in baselineDict.items()])/(1.0*precisionFactor))
np.savetxt(scriptDirectory + "/redundancy.dat",np.asarray([len(uniqueBaseline[1]) for uniqueBaseline in baselineDict.items()]),fmt='%i')

antennaPairDict = {}
for uniqueIndex in range(len(baselineDict.items())):
	allItems = baselineDict.items()
	for antennaPair in allItems[uniqueIndex][1]:
		antennaPairDict[antennaPair] = uniqueIndex
pickle.dump(antennaPairDict, open(scriptDirectory + "/antennaPairUniqueBaselineIndexDict.p", 'w'))