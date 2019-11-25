# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 17:01:15 2019

@author: bguilfoyle
"""

import numpy as np
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
from regions import PixCoord, EllipsePixelRegion
import matplotlib.pyplot as plt

x0, y0 = 15, 10
a, b = 8, 5
theta = Angle(30, 'deg')
e = Ellipse2D(amplitude=100., x_0=x0, y_0=y0, a=a, b=b, theta=theta.radian)
y, x = np.mgrid[0:20, 0:30]
fig, ax = plt.subplots(1, 1)
ax.imshow(e(x, y), origin='lower', interpolation='none', cmap='Greys_r')

center = PixCoord(x=x0, y=y0)
reg = EllipsePixelRegion(center=center, width=2*a, height=2*b, angle=theta)
patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
ax.add_patch(patch)