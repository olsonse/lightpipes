#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
Twyman Green interferometer.
The Twyman green interferometer consists of a beamsplitter
and two mirrors. One of the mirrors
is perfect, the other can be investigated for aberrations. 

F.A. van Goor, March 1998.
"""

import common
from lightpipes import Field
from pylab import *

m=1
cm=1e-2*m
mm=1e-3*m
nm=1e-9*m

wavelength  = 500*nm
size        = 30*mm
R           = 5*mm
N           = 250
z1          = 50*cm
z2          = 40*cm
z3          = 40*cm
z4          = 100*cm
RBS         = 0.3
nz          = 3
mz          = 1
Rz          = 0.005
Az          = 25

F1 = Field(N,size,wavelength)
F1.circular_aperture(R, 0, 0)
F1.forvard(z1)

F2 = F1.copy()

F1 *= RBS;          F2 *= 1-RBS
F1.forvard(2*z2);   F2.forvard(z3)
F1;                 F2.zernike(nz, mz, Rz, Az)
F1;                 F2.forvard(z3)
F1 *= RBS;          F2 *= 1-RBS

F = F1 + F2

F.forvard(z4)
F.interpolate(0.012, 250, 0, 0, 0, 1)

imshow( abs(F.value)**2 ) # plot intensity
title('Twyman Green interferometer with Zernike aberration')
show()
