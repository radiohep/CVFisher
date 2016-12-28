#! /usr/bin/env python

import sys, optparse, numpy as np
from scipy import interpolate

o = optparse.OptionParser()
o.set_usage('scalePspec.py [options]')
o.set_description(__doc__)
o.add_option('-R', '--restFreq', dest='restFreq', default=1.42040575177, type=float,
        help="Rest frequency of spectral line in question. In GHz")
o.add_option('-F', '--chosenFreq', dest='chosenFreq', default=0.800, type=float,
        help="Frequency in question. In GHz")
o.add_option('--pspecPath', dest='pspecPath', type=str,
        help="Path to power spectrum file")
o.add_option('--bgPath', dest='bgPath', type=str,
        help="Path to background file")
o.add_option('--outPath', dest='outPath', type=str,
        help="Where to put output Delta^2(k) file")
opts, args = o.parse_args(sys.argv[1:])

restFreq = opts.restFreq
chosenFreq = opts.chosenFreq
pspecPath = opts.pspecPath
bgPath = opts.bgPath
outPath = opts.outPath

chosenRedshift = restFreq / chosenFreq - 1


# Load matter power spectrum
matterPspec = np.loadtxt(pspecPath)
# Load file with global signal, bias etc.
bgData = np.loadtxt(bgPath)

temperature = interpolate.interp1d(bgData[:,0], bgData[:,3])(chosenRedshift)
growth = interpolate.interp1d(bgData[:,0], bgData[:,2])(chosenRedshift)
bias = interpolate.interp1d(bgData[:,0], bgData[:,4])(chosenRedshift)

Pz0 = interpolate.interp1d(matterPspec[:,0], matterPspec[:,1])
DeltaSqz = temperature**2 * growth**2 * bias**2 * matterPspec[:,1]
DeltaSqz *= matterPspec[:,0]**3 / ( 2. * np.pi**2 )

np.savetxt(outPath,np.array([matterPspec[:,0],DeltaSqz]).T)
