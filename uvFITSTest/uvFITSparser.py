import pyfits
import numpy as np

#UV_Fits_Data_Path='21cm_Poisson_PL_FD_PL_SCdist_0d063GN.uvfits'
UV_Fits_Data_Path='blind.uvfits'

###
#Open file and read in data
###
HDUlist=pyfits.open(UV_Fits_Data_Path, memmap=True)
Data=HDUlist[0].data
###

###
#Testing
###
print 'Nvis:', len(Data)
print 'Elements per Vis:', len(Data[0])
#print 'Nvis:', len(Data[0][-1])
print 'Per vis metadata:', len(Data[0][:-1])
print 'Nchan:', len(Data[0][-1][0][0][0])
print 'Polarizations:', len(Data[0][-1][0][0][0][0])
print 'len 1 pol real, imag, weight:', len(Data[0][-1][0][0][0][0][0])
###

print 'Generate array from fits file...'
Data=np.array(Data)
print 'Array generated. Continuing.'

#for key in HDUlist[0].header.keys():
#	print key, HDUlist[0].header[key]

#print Data[0]

#print HDUlist[0].header

#print HDUlist[0].header['TIMES']
#print Data[0][-1]


#print np.asarray([Data[y][6] for y in range(len(Data))])

print np.unique(np.asarray([Data[y][5] for y in range(len(Data))]))

visibilities = np.asarray([Data[y][9] for y in range(len(Data))])

visibilitiesForAnt258 = visibilities[np.asarray([Data[y][5] for y in range(len(Data))]) == 7197]
#baselineLength = 


#print len(visibilitiesForAnt258[0][0][0][0][0][0])
#print len(visibilitiesForAnt258[0][0][0][0][0])
#print len(visibilitiesForAnt258[0][0][0][0])
#print len(visibilitiesForAnt258[0][0][0])
#print len(visibilitiesForAnt258[0][0])
#print len(visibilitiesForAnt258[0])
visibilitiesForAnt258 = visibilitiesForAnt258[:,0,0,0,:,:,:]


import matplotlib.pyplot as plt

phasesForAnt258XX = np.arctan2(visibilitiesForAnt258[:,:,0,0],visibilitiesForAnt258[:,:,0,1])
print phasesForAnt258XX.shape
realsForAnt258XX = visibilitiesForAnt258[:,:,0,0]
plt.imshow(phasesForAnt258XX, interpolation='nearest')
plt.colorbar()
plt.show()




#indexing scheme = data[baselineNumber cycled with time][useless][useless][useless][freqChannel][polarization][real/imag/weight]

