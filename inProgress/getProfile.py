
"""
Created on Wed Feb 5 16:35:16 2020

@author: bguilfoyle - bengfoyle.github.io

@overview: Map n ellipses to an image, and return the brightness profile
"""

#astropy utilities for viewing and editing fits files
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import Angle
from astropy.io import fits
#matplot to map ellipses, and display images
import matplotlib.patches as mpatches
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
from matplotlib.path import Path
#Curve fitting for mapping Sersic profile
from scipy.optimize import curve_fit
#gui elements
from tkinter import *
#trig identidies (sin,cos)
import numpy as np

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
def sersic(r,I_e,b_n,r_e,n,offset):
    """
    Overview: Evaluate the Sersic profile for given parameters
    """
    return I_e * np.exp(-b_n*(pow((r/r_e),(1/n))-1)) + offset
#==============================================================================

#===============================================================================
def getProfile():
    """
    Overview: Construct ellipses around a user defined image, return the
    average brightness of each ellipse, and a brightness profile plot.
    """
    image_file = get_pkg_data_filename("finalCombined.fit")#file.get())
    image_data = fits.getdata(image_file, ext=0)
    mini = 0.005

    #Cretaing initial elipse parameters
    x0, y0 = 450,342#float(xCoord.get()),float(yCoord.get())
    a, b = 0.5,0.2#float(majorAxis.get()),float(minorAxis.get())
    theta = Angle(float(-50),"deg")#angleInput.get()), 'deg')

    y, x = np.mgrid[0:50, 0:50]
    fig, ax = plt.subplots(1, 1)
    #display galaxy image
    ax.imshow(loggray(image_data,mini), origin='lower', interpolation='none', \
    cmap='Greys_r')
    pathList = []
    major = []
    minor = []
    focus = []
    n = 25#int(numEllipse.get()) # number of ellipses
    numRegions = int(numEntry.get()) #number of regions to map
    growth = 5#float(growthFactor.get()) #factor by which each ellipse will grow

    """
    Map ellipses to image
    """
    for i in range(0,n):
        #matplot patch in the shape of ellipse
        shape = mpatches.Ellipse((x0, y0), 2*a, 2*b, theta.degree,\
                            edgecolor='red',facecolor='none')
        ax.add_patch(shape)
        pathList.append(list(shape.get_verts())) #points along edge of ellipse
        #Record dimensions of elipse, and increase by growth
        major.append(a)
        minor.append(b)
        focus.append(np.sqrt(a**2 - b**2)) #distance to ellipse focii
        a += a/growth; b += b/growth
    plt.show()


    """
    Iterate over each ellipse, for each point along that ellipse sum the
    """
    avgBrightness = []
    region = list(np.array_split(pathList,numRegions))

    for r in range(0,numRegions):
        prevAvg = 10 #prevAvg set to placeholder of
        for i in range(0,len(region[r])):
            e = region[r][i]
            bright = 0
            numPoints = len(e)
            for j in range(0,len(e)):
                #check for outliers by comparing to previous average brightness
                intensity = image_data[int(e[j][0])][int(e[j][1])]
                if intensity > prevAvg:
                    numPoints -= 1
                else:
                    bright += intensity
            try:
                avgBrightness.append(bright/numPoints)
                prevAvg = avgBrightness[i]
            except:
                pass

    avgBrightness = list(np.array_split(avgBrightness,numRegions))
    focus = list(np.array_split(focus,numRegions))
    b_n = bnEntry.get().split(",")
    b_n = [float(b) for b in b_n]
    r_e = reEntry.get().split(",")
    r_e = [float(r) for r in r_e]
    n = nEntry.get().split(",")
    n = [float(num) for num in n]

    for i in range(0,len(avgBrightness)):
        I_e = avgBrightness[i][0]
        x = np.linspace(0,focus[i][len(focus[i])-1],len(focus[i]))
        try:
            I_r = [sersic(xVal,I_e,b_n[i],r_e[i],n[i],0) for xVal in x]
            popt, pcov = curve_fit(sersic, focus[i], avgBrightness[i], p0 = (I_e,b_n[i],r_e[i],n[i],0))
            fittedData = sersic(x, *popt)
            name = "Region ",i
            print(i)
            plt.plot(focus[i],fittedData, label=name)
        except:
            i +=1

    #plot distance to focus as a function of brightness
    print("****")
    plt.plot(focus,avgBrightness,label = "Experimental Data")
    print("&&&&")
    #plt.plot(focus,I_r,label = "Sersic Profile")
    plt.legend()
    plt.show()
    print(popt)
#===============================================================================

#===============================================================================
"""
GUI Elements
"""
window = Tk()
window.title("Ellipse Mapping & Brightness Tool")
window.geometry('500x400')

#Insert text with user input text field.

filePrompt = Label(window, text="Filename.fit(s):")
filePrompt.grid(column=0, row=0)
file = Entry(window, width = 25)
file.grid(column = 1, row = 0)

origin = Label(window, text="Origin Point of Galaxy (x,y):")
origin.grid(column=0, row=1)

xCoord = Entry(window, width = 25)
xCoord.grid(column = 1, row = 1)
yCoord = Entry(window, width = 25)
yCoord.grid(column = 2, row = 1)

minmajPrompt = Label(window, text="Semimajor/minor axis of ellipse:")
minmajPrompt.grid(column=0, row=2)

majorAxis = Entry(window, width = 25)
majorAxis.grid(column = 1, row = 2)
minorAxis = Entry(window, width = 25)
minorAxis.grid(column = 2, row = 2)

numEllipsePrompt = Label(window, text="Number of ellipses:")
numEllipsePrompt.grid(column = 0, row = 3)
numEllipse = Entry(window, width = 25)
numEllipse.grid(column = 1, row = 3)

growthPrompt = Label(window, text="Factor of growth:")
growthPrompt.grid(column = 0, row = 4)
growthFactor = Entry(window, width = 25)
growthFactor.grid(column = 1, row = 4)

anglePrompt = Label(window, text = "Angle in degrees")
anglePrompt.grid(column = 0, row = 5)
angleInput = Entry(window, width = 25)
angleInput.grid(column = 1, row = 5)

"""
Sersic parameters
"""

bnLabel = Label(window, text = "bn parameter")
reLabel = Label(window, text = "re parameter")
nLabel = Label(window, text = "n parameter")
numRegions = Label(window, text = "Number of Regions")
bnEntry = Entry(window, width = 25)
reEntry = Entry(window, width = 25)
nEntry = Entry(window, width = 25)
numEntry = Entry(window, width = 25)

bnLabel.grid(column = 0, row = 7)
reLabel.grid(column = 0, row = 8)
nLabel.grid(column = 0, row = 9)
numRegions.grid(column = 0, row = 10)
bnEntry.grid(column = 1, row = 7)
reEntry.grid(column = 1, row = 8)
nEntry.grid(column = 1, row = 9)
numEntry.grid(column = 1, row = 10)
submit = Button(window, text="Run with these Parameters", command = getProfile)
submit.grid(column=0, row=11)

#loop until closed
window.mainloop()
#===============================================================================
