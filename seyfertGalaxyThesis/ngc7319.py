# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:05:05 2019

@author: bguilfoyle
"""
import numpy as np
import astropy.io.fits as pyfits
from astropy.nddata import Cutout2D
import glob
import matplotlib.pyplot as plt
from photutils.isophote import EllipseGeometry
from photutils import EllipticalAperture
from photutils.isophote import Ellipse
from photutils.isophote import build_ellipse_model

#==============================================================================
def avgImage(img):
    """
    Overview: return average of an image
    """
    imageConcat = []
    for image in img:
        imageConcat.append(pyfits.getdata(image))

    finalImage = np.zeros(shape=imageConcat[0].shape)

    for image in imageConcat:
        finalImage += image
    return finalImage
#==============================================================================


#==============================================================================
def loggray(x, a=None, b=None):
    """
    Overview: Auxiliary function that specifies the logarithmic gray scale.
    a and b are the cutoffs : if not specified, min and max are used
    """
    if a == None:
        a = np.min(x)
    if b == None:
        b = np.max(x)
    linval = 10.0 + 990.0 * (x-float(a))/(b-a)
    return (np.log10(linval)-1.0)*0.5 * 255.0
#==============================================================================
path = "C:/Users/bguilfoyle/Documents/CompPhysics/FYP/seyfertGalaxyThesis/finalR.fit"
img = avgImage(glob.glob(path))
#plt.imshow(img)
position = (450,343)
size = (30, 50)
geometry = EllipseGeometry(x0=450, y0=343, sma=20, eps=0.5, pa=20.*np.pi/180.)
aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,geometry.\
                          sma*(1 - geometry.eps),geometry.pa)
ellipse = Ellipse(img, geometry)
isolist = ellipse.fit_image()
model_image = build_ellipse_model(img.shape, isolist)
residual = img - model_image

cutout = Cutout2D(img, position, size)

#hdu = pyfits.PrimaryHDU(cutout)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("NGC7319Cutout.fit")
#hdulist.close()



plt.imshow(cutout, origin='lower')
aper.plot(color='white')

#plt.imshow(cutout.data)
#plt.colorbar()