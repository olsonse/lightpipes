#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
March 1998. F.A. van Goor.
Two Holes Interferometer.
"""

import common
from lightpipes import Field
from pylab import *

m  = 1
nm = 1e-9*m
mm = 1e-3*m
cm = 1e-2*m

wavelength  = 550*nm
size        = 5*mm
N           = 300
R           = 0.12*mm
d           = 0.5*mm
z           = 5*cm

F = Field(N,size,wavelength)
F2 = F.copy()
F.circular_aperture(R,d,0)
F2.circular_aperture(R,-d,0)
F.value[:] = F.value + F2.value
del F2

for l in xrange(1,11):
  F.forvard(z)
  subplot(2,5,l)
  imshow( abs(F.value)**2 )
  xticks([]), yticks([])
  title( 'z={z} cm'.format(z =(l*z/cm) ) )

show()
