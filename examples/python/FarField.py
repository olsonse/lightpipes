#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
Calculates the Fraunhofer diffraction of a round hole.
"""

import common
from lightpipes import Field
from pylab import *

m=1
nm=1e-9*m
mm=1e-3*m
cm=1e-2*m

wavelength  = 1000*nm
size        = 10*cm
N           = 150
R           = 10*mm
z           = 1000*m
f1          = 200*m
f2          = -200*m

F = Field(N,size,wavelength)
F.circular_aperture(R)
F.lens(f1)
F.lens_forvard(f2, z)

imshow( abs(F.value)**2 ) # plot intensity
show()
