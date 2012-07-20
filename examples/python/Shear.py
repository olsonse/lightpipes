#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
Shearing interferometer with aberrated wavefront.
"""

import common
from lightpipes import Field
from pylab import *

m=1
nm=1e-9*m
mm=1e-3*m
cm=1e-2*m

wavelength  = 500*nm
size        = 4*cm
N           = 256
R           = 1*cm
f           = -20*m
z1          = 50*cm
z2          = 50*cm
D           = 3*mm
D1          = 1*mm
Rplate      = 0.5

nZ = [4,3,2,0]
mZ = [0,-1,2,0]
RZ = [10*mm,10*mm,7*mm,0]
AZ = [10,10,10,0]

choices = [
  "Spherical aberration with: n=4, m=0, R=10mm, A=10",
  "Coma with: n=3, m=-1, R=10mm, A=10",
  "Astigmatism with: n=2, m=2, R=7mm, A=10",
  "No aberration",
]

while True:
  print """
    CHOOSE ABERRATION
  0.  {0}
  1.  {1}
  2.  {2}
  3.  {3}
  <other>  STOP
  """.format(*choices)
  choice = raw_input()
  if choice not in ['0','1','2','3']:
    break
  choice = int(choice)


  F = Field(N,size,wavelength)
  F.circular_aperture(R)
  if choice == 3:
    F.lens(f)
  else:
    k = choice
    F.zernike(nZ[k],mZ[k],RZ[k],AZ[k])

  F.forvard(z1)
  F1 = F.copy()
  F1.value[:] *= Rplate
  F.value[:] *= (1 - Rplate)

  F.interpolate(size,N,D,D1)
  F.value[:] += F1.value
  Int = (F.value * F.value.conj()).real

  imshow(Int)
  title( choices[choice] )
  show()
