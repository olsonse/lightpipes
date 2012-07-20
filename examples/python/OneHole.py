#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
One Hole Diffraction.
"""
import common
from lightpipes import Field
from pylab import *

m=1
nm=1e-9*m
mm=1e-3*m
cm=1e-2*m

wavelength  = 550*nm;
size        = 5*mm
N           = 100
R           = 1*mm
z           = 25*cm

F = Field(N,size,wavelength)
F.circular_aperture(R)
F.rectangular_screen(R*2, R/8, 0, 0, -45)
F.circular_screen(R/4)
F.fresnel(z)

imshow( abs(F.value)**2 ) # plot intensity
title('Intensity Distribution in the somewhat far-field');
show()
