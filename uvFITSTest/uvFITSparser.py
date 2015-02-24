import pyfits
import numpy as np

UV_Fits_Data_Path='21cm_Poisson_PL_FD_PL_SCdist_0d063GN.uvfits'

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

print [Data[y][6] for y in range(10000,12000,1)]
