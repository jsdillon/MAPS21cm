import numpy as np
import healpy as hp
import math
import os
import matplotlib
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
import ephem
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
from Source.GlobalSkyModel import GlobalSkyModel
import scipy.constants as const
import cPickle as pickle
import os
from Mapmaker import Mapmaker

plt.close("all")
resultsFolder = Mapmaker(overridePSFExtension = 5)

#Load all results
s = pickle.load(open(resultsFolder + "specifications.p","rb"))
coords = pickle.load(open(s.resultsFolder + "coordinates.p","rb"))
times = pickle.load( open(s.resultsFolder + "times.p","rb"))
ps = pickle.load(open(s.resultsFolder + "pointSourceCatalog.p","rb"))
Dmatrix = np.load(s.resultsFolder + "Dmatrix.npy")
PSF = np.load(s.resultsFolder + "PSF.npy")
coaddedMap = np.load(s.resultsFolder + "coaddedMap.npy")
mapNoiseCovariance = np.load(s.resultsFolder + "mapNoiseCovariance.npy")
if s.PSFforPointSources and ps.nSources>0:
    pointSourcePSF = np.load(s.resultsFolder + "pointSourcePSF.npy")


s.GSMNSIDE = s.mapNSIDE
coordsGSM = Geometry.Coordinates(s,True)
GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))
convolvedGSM = np.dot(PSF,interpoltedGSMRotated[coords.extendedIndices])

def plotFacet(s,coords,facetMap,plotTitle):
    mapToPlot = np.zeros(coords.mapPixels)
    mapToPlot[coords.mapIndices] = facetMap    
    hp.mollview(mapToPlot, title=plotTitle)

plotFacet(s,coords,convolvedGSM,"Convolved GSM")
plotFacet(s,coords,coaddedMap,"coaddedMap")

#hp.mollview(interpoltedGSMRotated, title='test')



