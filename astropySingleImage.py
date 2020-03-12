# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:05:05 2019

@author: bguilfoyle
"""
import numpy as np
import astropy.io.fits as pyfits
import glob
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from tkinter import *

#==============================================================================
def lingray(x, a=None, b=None):
    """
    Overview: Auxiliary function that specifies the linear gray scale.
    a and b are the cutoffs : if not specified, min and max are used
    """
    if a == None:
        a = np.min(x)
    if b == None:
        b = np.max(x)
    return 255.0 * (x-float(a))/(b-a)
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

#==============================================================================
def correctedImage(raw,dark,flat,bias):
    """
    Overview: Return a corrected image based off basic raw reduction
    """
    return (raw - bias - dark) / flat
#==============================================================================

def addPlot(image,colour,newAlpha):
    """
    Overview: Make a plot using plt.imshow
    """
    plt.imshow(loggray(image), cmap= colour, alpha = newAlpha)
    plt.colorbar()
#==============================================================================    

path = "C:/Users/bguilfoyle/Documents/CompPhysics/FYP/seyfertGalaxyThesis/data/seyfertImages/"
rawRPath = path + "raw/*R.fit"
rawHPath = path + "raw/*H.fit"
rawSIIPath = path + "raw/*SII.fit"
rawVPath = path + "raw/*V.fit"
rawPath = path + "raw/*.fit"
biasPath = path + "bias/*.fit"
darkPath = path + "dark/*bin2.fit"
flatPath = path + "flat/Flat_bin2*.fit"

raw = glob.glob(rawRPath)
bias = glob.glob(biasPath)
dark = glob.glob(darkPath)
flat = glob.glob(flatPath)

avgRaw = avgImage(raw)
avgBias = avgImage(bias)
avgDark = avgImage(dark)
avgFlat = avgImage(flat)

finalCombined = correctedImage(avgRaw,avgDark,avgFlat,avgBias)
plt.imshow(loggray(finalCombined), cmap = "Reds")
plt.colorbar()