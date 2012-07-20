#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
calculation of the intensity near the
focus of a lens with lightpipes.Field.steps.
F.A. van Goor
"""

import common
from lightpipes import Field
from pylab import *

m  = 1
nm = 1e-9*m
mm = 1e-3*m
cm = 1e-2*m

wavelength  = 632.8*nm
size        = 4*mm
N           = 100
R           = 1.5*mm
dz          = 10*mm
f           = 50*cm
n           = (1 + .1j)*ones( (N,N) )

F = Field(N,size,wavelength)
F.circular_aperture(R,0,0)
F.lens(f,0,0)

Icross = zeros( (N,N) )
for l in xrange(N):
  F.steps(dz,1,n)
  Int = (F.value * F.value.conj()).real

  for k in xrange(N):
    Icross[l,k] = Int[N/2,k]
 

imshow(Icross)
show()
