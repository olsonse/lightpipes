#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
"""

from common import Propagate
from lightpipes import Field, lg
from pylab import *

m=1
nm=1e-9*m
um=1e-6*m
mm=1e-3*m
cm=1e-2*m

N               = 512                   # Number of pixels
side_length     = N* 15*um              # physical size of SLM
wavelength      = 780*nm                # Wavelength
gaussian_size   = 1.00/sqrt(2)*mm       # 1/e Intensity width of beam; (note, it's not the 1/e^2)
axicon_angle    = 175/180.*pi           # included angle of SLM
axicon_n1       = 1.5                   # index of refraction of the axicon
step_size       = 5.*mm                 # LPForvard step size
l               = 2                     # azimuthal order of LG beam

######## END Param ##############


####### Create LG order 
dx = (2/(N-1.))
x = y = 3*r_[-1:(1+dx):dx]
[xx, yy] = meshgrid(x,y)

lg = lg.LG_xy( l=l, p=0, xx=xx, yy=yy, omega0=2/sqrt(l) )
# normalize lg; we only want the phase information
lg /= abs(lg)
#imshow (x,y,arg(lg))
####### END Create LG order 



propagate = Propagate()

F = Field( N, side_length, wavelength )
F.gaussian_aperture(gaussian_size)
F.axicon(axicon_angle, axicon_n1)
#propagate(F, z=13.*cm, dz=step_size)
F.forvard(7.5*cm)
F.axicon(axicon_angle, axicon_n1)

print 'moving towards slm'
#F.forvard(40*cm)
propagate(F, z=75.*cm, dz=5*cm)

print 'first lens'
F.lens(-75*cm)
propagate(F, z=100.*cm, dz=5*cm)

print 'second lens'
F.lens(-25*cm)
propagate(F, z=25.*cm, dz=1*cm)

print 'applying slm'
F.value[:] *= lg
propagate(F, z=50.*cm, dz=2*cm)

print 'press enter to exit'
raw_input()
