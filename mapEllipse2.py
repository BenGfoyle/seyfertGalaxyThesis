
"""
Created on Wed Feb 5 16:35:16 2020

@author: bguilfoyle - bengfoyle.github.io

@overview: Map n ellipses to an image using astropy.
"""


from astropy.utils.data import get_pkg_data_filename
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
from astropy.io import fits

import matplotlib.patches as mpatches
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
from matplotlib.path import Path
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

image_file = get_pkg_data_filename('finalCombined.fit')
image_data = fits.getdata(image_file, ext=0)
print(image_data.shape)
mini = 0.005

#Cretaing initial elipse parameters
x0, y0 = 450, 340
a, b = 25, 10
theta = Angle(-50, 'deg')

y, x = np.mgrid[0:50, 0:50]
fig, ax = plt.subplots(1, 1)
#display galaxy image
ax.imshow(loggray(image_data,mini), origin='lower', interpolation='none', \
cmap='Greys_r')

n = 10 # number of ellipses
#plot number of elipses at various distances appart
for i in range(0,n):
    shape = mpatches.Ellipse((x0, y0), 2*a, 2*b, theta.degree, edgecolor='red',
                          facecolor='none')
    ax.add_patch(shape)
    elpath=list(shape.get_verts())
    print(elpath)
    #Increase size of elipse for next iteration
    a += a/5; b += b/5

plt.show()
