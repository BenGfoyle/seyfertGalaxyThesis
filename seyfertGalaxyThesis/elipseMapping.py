# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:00:54 2019

@author: bguilfoyle
"""

from photutils.isophote import EllipseGeometry
import numpy as np
from imageReduction import loggray
import glob
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

geometry = EllipseGeometry(x0=75, y0=75, sma=20, eps=0.5, pa=20.*np.pi/180.)

plt.imshow(loggray(glob.glob("finalR.fit")))
plt.colorbar()