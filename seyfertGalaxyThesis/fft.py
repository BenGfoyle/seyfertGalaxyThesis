# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:30:43 2019

@author: bguilfoyle
"""

from scipy import fftpack
import astropy.io.fits as pyfits
import numpy as np
import pylab as py
import radialProfile
image = pyfits.getdata("finalCombined.fit")

# Take the fourier transform of the image.
F1 = fftpack.fft2(image)

# Now shift the quadrants around so that low spatial frequencies are in
# the center of the 2D fourier transformed image.
F2 = fftpack.fftshift( F1 )

# Calculate a 2D power spectrum
psd2D = np.abs( F2 )**2

# Calculate the azimuthally averaged 1D power spectrum
psd1D = radialProfile.azimuthalAverage(psd2D)

# Now plot up both
py.figure(1)
py.clf()
py.imshow( np.log10( image ), cmap=py.cm.Greys)

py.figure(2)
py.clf()
py.imshow( np.log10( psd2D ))

py.figure(3)
py.clf()
py.semilogy( psd1D )
py.xlabel("Spatial Frequency")
py.ylabel("Power Spectrum")

py.show()