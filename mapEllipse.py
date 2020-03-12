"""
Created on Wed Feb 5 16:35:16 2020

@author: bguilfoyle - bengfoyle.github.io

@overview: Map n ellipses to an image using astropy, and return the average
brightness along each ellipse
"""


from astropy.utils.data import get_pkg_data_filename
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
from astropy.io import fits

import matplotlib.patches as mpatches
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import scipy.optimize

import pandas as pd

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
def anglesInEllipse(num,a,b):
    assert(num > 0)
    assert(a < b)
    angles = (2 * np.pi * np.arange(num) / num)
    if a != b:
        e = (1.0 - a ** 2.0 / b ** 2.0) ** 0.5
        tot_size = sp.special.ellipeinc(2.0 * np.pi, e)
        arc_size = tot_size / num
        arcs = np.arange(num) * arc_size
        res = sp.optimize.root(
            lambda x: (sp.special.ellipeinc(x, e) - arcs), angles)
        angles = [(x + 0.87)% 2*np.pi for x in angles]
        print(angles)
    return angles
#==============================================================================

#===============================================================================
def radialDistance(x,y):
    x0, y0 = 450, 340
    return np.sqrt((x - x0)**2 + (y - y0)**2)
#===============================================================================

image_file = get_pkg_data_filename('finalCombined.fit')
image_data = fits.getdata(image_file, ext=0)
mini = 0.005

#Cretaing initial elipse parameters
x0, y0 = 450, 340
a, b = 1, 2.5
theta = 1.57

fig, ax = plt.subplots(1, 1)
#display galaxy image
image = loggray(image_data,mini)
ax.imshow(image, origin='lower', interpolation='none', cmap='Greys_r')

n = 15 #number of ellipses
nPoints = 100
ePoints = [] #store the points present in each ellipse
avgBrightnes = []
#plot number of elipses at various distances appart
for i in range(0,n):
    phi = anglesInEllipse(nPoints,a,b)
    e = (1.0 - a ** 2.0 / b ** 2.0) ** 0.5
    arcs = sp.special.ellipeinc(phi, e)
    xAxis = x0 + (b * np.sin(phi)); yAxis = y0 + (a * np.cos(phi))
    ePoints.append([xAxis,yAxis])
    ax.plot(xAxis,yAxis)
    print(xAxis,yAxis,phi)

    brightness = [image[int(xAxis[i])][int(yAxis[i])] for i in range(0,len(xAxis))]
    avgBrightnes.append(sum(brightness) / nPoints)

    a += a/3; b += b/3
#print(avgBrightnes)

# data = pd.DataFrame(ePoints)
# data.to_csv('ellipsePoints.csv')
plt.show()
