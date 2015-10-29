#!/usr/bin/env python

import aipy as a, numpy as n, pylab as p, pyfits
import sys, os, optparse
from mpl_toolkits.basemap import Basemap

class HealpixMap(a.healpix.HealpixMap):
    '''Subclass to match the slightly different Healpix format James wrote.'''
    def from_fits(self, filename, hdunum=1, colnum=0):
        hdu = pyfits.open(filename)[hdunum]
        data = hdu.data.field(colnum)
        data = data.flatten()
        if not data.dtype.isnative:
            data.dtype = data.dtype.newbyteorder()
            data.byteswap(True)
        scheme= hdu.header['ORDERING'][:4]
        self.set_map(data, scheme=scheme)

def data_mode(data, mode='abs'):
    '''Convert data to desired plotting mode.'''
    if mode.startswith('phs'): data = n.angle(data)
    elif mode.startswith('lin'):
        data = n.absolute(data)
        data = n.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = n.absolute(data)
        data = n.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data


o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True,cal=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('-i', '--interp',  action='store_true', 
    help="Interpolate")

opts,args = o.parse_args(sys.argv[1:])
cmap = p.get_cmap(opts.cmap)

map = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)

for filename in args:
    print 'Reading', filename
    h = HealpixMap(fromfits= [0])
    print 'SCHEME:', h.scheme()
    print 'NSIDE:', h.nside()
    if opts.interp: h.set_interpol(True)

    lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
    lons = 360 - lons
    lats *= a.img.deg2rad; lons *= a.img.deg2rad
    y,x,z = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
    try: data, indices = h[x,y,z]
    except(ValueError): data = h[x,y,z]
    data.shape = lats.shape
    data /= h[0,0,1]
    #data = data**2 # only if a voltage beam
    data = data_mode(data, opts.mode)

    map.drawmapboundary()
    map.drawmeridians(n.arange(0, 360, 30))
    map.drawparallels(n.arange(0, 90, 10))

    if opts.max is None: max = data.max()
    else: max = opts.max
    if opts.drng is None: min = data.min()
    else: min = max - opts.drng
    step = (max - min) / 10
    levels = n.arange(min-step, max+step, step)
    map.imshow(data, vmax=max, vmin=min, cmap=cmap)
    #map.contourf(cx,cy,data,levels,linewidth=0,cmap=cmap)
    p.colorbar()

    p.show()
