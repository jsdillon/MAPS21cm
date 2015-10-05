import numpy as np
import healpy as hp
import math
import os
import ephem
import sys
import matplotlib
matplotlib.use("pdf") 
import matplotlib.pyplot as plt
from Source.Specifications import Specifications
from Source.PrimaryBeams import PrimaryBeams
from Source.VisibilitySimulator import VisibilitySimulator
from Source import Geometry
from Source.PointSourceCatalog import PointSourceCatalog
from Source import MatricesForMapmaking as MapMats
from Source.LoadVisibilities import LoadVisibilities
import scipy.constants as const
from Source.GlobalSkyModel import GlobalSkyModel
from scipy import interp
plt.close('all')

freq = 150
#Load in everything we need, figure out which LSTs to work with
print "Now working on calculating the observing season at " + str(freq) + " MHz..."             
mainDirectory = os.path.dirname(os.path.abspath(__file__))    

configFile = "configuration.txt"
s = Specifications(mainDirectory, configFile, freq)
s.simulateVisibilitiesWithPointSources = False
s.OverrideSpecifications({})

times = Geometry.Times(s)
#times.CutOutUnusedLSTsAndGroupIntoSnapshots(s)
coords = Geometry.Coordinates(s)
PBs = PrimaryBeams(s)
ps = PointSourceCatalog(s,times)

#visibilities = VisibilitySimulator(s,PBs,ps,times,coords)
#plt.semilogy(times.LSTs, np.abs(visibilities))


beamWeightedPower = np.zeros(len(times.LSTs))
if s.GSMNSIDE < s.mapNSIDE:
    s.GSMNSIDE = s.mapNSIDE
coordsGSM = Geometry.Coordinates(s,True)
GSM = GlobalSkyModel(s.freq, s.GSMlocation, s.GSMNSIDE)
interpoltedGSMRotated = hp.get_interp_val(GSM.hpMap,-coordsGSM.galCoords.b.radian+np.pi/2, np.asarray(coordsGSM.galCoords.l.radian))        

#loop over times and baselines to calculate visibilities
for t in range(len(times.LSTs)):
    pixelAlts, pixelAzs = Geometry.convertEquatorialToHorizontal(s,coordsGSM.pixelRAs,coordsGSM.pixelDecs,times.LSTs[t])
    rHatVectors = Geometry.convertAltAzToCartesian(pixelAlts,pixelAzs)
    primaryBeam = hp.get_interp_val(PBs.beamSquared("X","x",s.pointings[t]), np.pi/2-pixelAlts, pixelAzs)
    beamWeightedPower[t] = np.sum(interpoltedGSMRotated * primaryBeam)  / (np.sum(primaryBeam))

#    plt.semilogy(times.LSTs, beamWeightedPower)


HSA = ephem.Observer()
HSA.lat = str(s.arrayLat)
HSA.lon = str(s.arrayLong)
HSA.date = '2016/1/1 12:0:0'
HSA2 = ephem.Observer()
HSA2.lat = str(s.arrayLat)
HSA2.lon = str(s.arrayLong)
sun = ephem.Sun()
sunrises = []
sunsets = []
midnights = []
midnightLSTs = []
midnightDays = []
sunsetLSTs = []
sunsetDays = []
sunriseLSTs = []
sunriseDays = []

sunriseUTCs = []
midnightUTCs = []
sunsetUTCs = []
for day in range(365):
    HSA.date += 1
    sunrise = HSA.next_rising(sun).tuple()
    sunrises.append(sunrise[3] + sunrise[4]/60.0 + sunrise[5]/3600.0)
    sunset = HSA.next_setting(sun).tuple()
    sunsets.append(sunset[3] + sunset[4]/60.0 + sunset[5]/3600.0)
    midnight = ephem.Date(HSA.next_setting(sun) + (HSA.next_rising(sun) - HSA.next_setting(sun))/2).tuple()
    midnights.append(midnight[3] + midnight[4]/60.0 + midnight[5]/3600.0)
    HSA2.date = ephem.Date(HSA.next_setting(sun) + (HSA.next_rising(sun) - HSA.next_setting(sun))/2)
    midnightLSTs.append(float(HSA2.sidereal_time()) * 12.0 / np.pi)
    midnightDays.append(ephem.Date(HSA.next_setting(sun) + (HSA.next_rising(sun) - HSA.next_setting(sun))/2) - ephem.Date('2016/1/1 0:0:0'))
    
    HSA2.date = HSA.next_setting(sun)
    sunsetLSTs.append(float(HSA2.sidereal_time()) * 12.0 / np.pi)
    sunsetDays.append(ephem.Date(HSA.next_setting(sun)) - ephem.Date('2016/1/1 0:0:0'))
    
    HSA2.date = HSA.next_rising(sun)
    sunriseLSTs.append(float(HSA2.sidereal_time()) * 12.0 / np.pi)
    sunriseDays.append(ephem.Date(HSA.next_setting(sun)) - ephem.Date('2016/1/1 0:0:0'))

       
    

#    plt.plot(np.asarray(sunrises))
#    plt.plot(np.asarray(sunsets))
#    allMidnights = np.transpose(np.asarray([midnightDays, midnightLSTs]))
#    allMidnights = allMidnights[allMidnights[:,1].argsort()]
#    daysMidnight = interp1d(allMidnights[:,1], allMidnights[:,0], kind='linear', bounds_error= False, fill_value = -1)(times.LSTs)
#    beamWeightedPowerMidnight = beamWeightedPower[daysMidnight>0]    
#    daysMidnight = daysMidnight[daysMidnight>0]
#    
#    allSunsets = np.transpose(np.asarray([sunsetDays, sunsetLSTs]))
#    allSunsets = allSunsets[allSunsets[:,1].argsort()]    
#    daysSunset = interp1d(allSunsets[:,1], allSunsets[:,0], kind='linear', bounds_error= False, fill_value = -1)(times.LSTs)
#    beamWeightedPowerSunset = beamWeightedPower[daysSunset>0]    
#    daysSunset = daysSunset[daysSunset>0]
#    
#
#    allSunrises = np.transpose(np.asarray([sunriseDays, sunriseLSTs]))
#    allSunrises = allSunrises[allSunrises[:,1].argsort()]
#    daysSunrise = interp1d(allSunrises[:,1], allSunrises[:,0], kind='linear', bounds_error= False, fill_value = -1)(times.LSTs)
#    beamWeightedPowerSunrise = beamWeightedPower[daysSunrise>0]    
#    daysSunrise = daysSunrise[daysSunrise>0]
#    
#    
#
#    sunset, = plt.plot(daysSunset/30, beamWeightedPowerSunset,'o',label='Sunset')    
#    midnight, = plt.plot(daysMidnight/30, beamWeightedPowerMidnight,'o',label='Midnight')
#    sunrise, = plt.plot(daysSunrise/30, beamWeightedPowerSunrise,'o',label='Sunrise')
#
#    plt.legend(handles=[sunset, midnight, sunrise])
#    plt.xlabel('Months since January 1, 2016')
#    plt.ylabel('Beam Weighted Power at Sunset/Midnight/Sunrise')
#    plt.title('Observing Seasons for 10 Degree FWHM Gaussian Beam')
#    
HSA = ephem.Observer()
HSA.lat = str(s.arrayLat)
HSA.lon = str(s.arrayLong)
HSA.date = '2016/1/1 0:0:0'
LSTs = np.zeros((365,60*24))
for day in range(365):
    for minute in range(60*24):
        HSA.date += ephem.minute
        LSTs[day,minute] = float(HSA.sidereal_time()) * 12.0 / np.pi

powers = np.zeros((365,60*24))
for day in range(365):
    powers[day,:] = interp(LSTs[day,:], times.LSTs, beamWeightedPower)

#%% 

fig, ax = plt.subplots()
cax = ax.pcolorfast(np.arange(365)*12.0/365.0, np.arange(60*24)/60.0,np.transpose(powers),vmin = 0)
#, aspect='auto', interpolation = 'none', origin = 'lower',vmin = 0)
plt.xlabel('Months since January 1, 2016')
plt.ylabel('Time UTC (hours)')
sunset, = plt.plot(np.arange(365)*12.0/365.0, sunsets,'r-',label='Sunset',lw=3)    
midnight, = plt.plot(np.arange(365)*12.0/365.0, midnights,'y-',label='Midnight',lw=3)    
sunrise, = plt.plot(np.arange(365)*12.0/365.0, sunrises,'g-',label='Sunrise',lw=3)    
cbar = fig.colorbar(cax)
cbar.set_label('Sum of Beam-Weighted GSM in (K) / Sum of Beam')
ax.legend(handles=[sunset, midnight, sunrise])
ax.set_ylim([0, 24])
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)

#cbar.ax.set_ylim([0, np.amax(powers)])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar

plt.title('Observing Season for 10$^\circ$ FWHM Gaussian Beam')
plt.savefig('observing_seasons.pdf',format='pdf')



