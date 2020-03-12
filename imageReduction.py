# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:05:06 2019

@author: bguilfoyle - bengfoyle.github.io
"""
import numpy as np
import astropy.io.fits as pyfits
import glob
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.modeling.models import Gaussian2D
from tkinter import *

#==============================================================================

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
def addPlot(image):
    """
    Overview: Make a plot using plt.imshow
    """
    plt.imshow(loggray(image))
    plt.colorbar()
#==============================================================================

#==============================================================================
def correctedImage(raw,dark,flat,bias):
    """
    Overview: Return a corrected image based off basic raw reduction
    """
    return (raw - bias - dark) / flat
#==============================================================================

#==============================================================================
"""
Overview: File paths and averaging images
"""
path = "C:/Users/bguilfoyle/Documents/CompPhysics/FYP/seyfertGalaxyThesis/data/seyfertImages/"
rawRPath = path + "raw/*R.fit"
rawHPath = path + "raw/*H.fit"
rawSIIPath = path + "raw/*SII.fit"
rawVPath = path + "raw/*V.fit"
rawPath = path + "raw/*.fit"
biasPath = path + "bias/*.fit"
darkPath = path + "dark/*bin2.fit"
flatPath = path + "flat/Flat_bin2*.fit"

raw = glob.glob(rawPath)
bias = glob.glob(biasPath)
dark = glob.glob(darkPath)
flat = glob.glob(flatPath)

avgRaw = avgImage(raw)
avgBias = avgImage(bias)
avgDark = avgImage(dark)
avgFlat = avgImage(flat)

rawR = glob.glob(rawRPath)
rawH = glob.glob(rawHPath)
rawSII = glob.glob(rawSIIPath)
rawV = glob.glob(rawVPath)

avgRawR = avgImage(rawR)
avgRawH = avgImage(rawH)
avgRawSII = avgImage(rawSII)
avgRawV = avgImage(rawV)

finalR = correctedImage(avgRawR,avgDark,avgFlat,avgBias)
finalH = correctedImage(avgRawH,avgDark,avgFlat,avgBias)
finalSII = correctedImage(avgRawSII,avgDark,avgFlat,avgBias)
finalV = correctedImage(avgRawV,avgDark,avgFlat,avgBias)
finalCombined = correctedImage(avgRaw,avgDark,avgFlat,avgBias)

plt.imshow(finalV)

#finalR.writeto('finalR.fit')
#finalH.writeto('finalH.fit')
#finalSII.writeto('finalSII.fit')
#finalV.writeto('finalV.fit')
#finalCombined.writeto('finalCombined.fit')


#hdu = pyfits.PrimaryHDU(finalR)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("finalR.fit")
#hdu = pyfits.PrimaryHDU(finalH)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("finalH.fit")
#hdu = pyfits.PrimaryHDU(finalSII)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("finalSII.fit")
#hdu = pyfits.PrimaryHDU(finalV)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("finalV.fit")
#hdu = pyfits.PrimaryHDU(finalCombined)
#hdulist = pyfits.HDUList([hdu])
#hdulist.writeto("finalCombined.fit")
#hdulist.close()

#==============================================================================
#"""
#GUI for easier tuning
#"""
#def clicked(x = None): #perform opertion on button click
#    image = [finalR, finalH, finalSII, finalV, finalCombined]
#    alpha = [txt1.get(),txt2.get(),txt3.get(),txt4.get(),txt5.get()]
#    for i in range(0,len(alpha)):
#        if not alpha[i]:
#            alpha[i] = 0
#
#    alpha = [float(i) for i in alpha]
#    plt.close()
#
#    colour = ["Reds","Blues","Greens","Purples","Greys"]
#    for i in range(0,len(alpha)):
#        addPlot(image[i],colour[i],alpha[i])
#
##Define window, name, and parameters
#window = Tk()
#window.title("Image Tuner")
#window.geometry('450x300')
#
##Insert text with user input text field.
#lbl1 = Label(window, text="Alpha of R Filter")
#lbl1.grid(column=0, row=0)
#
#txt1 = Scale(window,from_ = 0, to_ = 1 ,width = 10, orient = HORIZONTAL, resolution = 0.01)
#txt1.grid(column = 1, row = 0)
#
#lbl2 = Label(window, text="Alpha of H Filter")
#lbl2.grid(column=0, row=1)
#
#txt2 = Scale(window,from_ = 0, to_ = 1 ,width = 10, orient = HORIZONTAL, resolution = 0.01)
#txt2.grid(column = 1, row = 1)
#
#lbl3 = Label(window, text="Alpha of S-II Filter")
#lbl3.grid(column=0, row=2)
#
#txt3 = Scale(window,from_ = 0, to_ = 1 ,width = 10, orient = HORIZONTAL, resolution = 0.01)
#txt3.grid(column = 1, row = 2)
#
#lbl4 = Label(window, text="Alpha of V Filter")
#lbl4.grid(column=0, row=3)
#
#txt4 = Scale(window,from_ = 0, to_ = 1 ,width = 10, orient = HORIZONTAL, resolution = 0.01)
#txt4.grid(column = 1, row = 3)
#
#lbl5 = Label(window, text="Alpha of Combined Filter Image")
#lbl5.grid(column=0, row=4)
#
#txt5 = Scale(window,from_ = 0, to_ = 1 ,width = 10, orient = HORIZONTAL, resolution = 0.01)
#txt5.grid(column = 1, row = 4)
#
##button that runs "clicked" function on click
#btn1 = Button(window, text="Submit Values", command = clicked)
#btn1.grid(column=1, row=5)
#
##loop until closed
#window.mainloop()
#==============================================================================
