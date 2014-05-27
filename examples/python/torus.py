#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox

*** Not sure if I completed this example...
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
side_length     = N* 15.*um             # physical size of SLM
wavelength      = 780*nm                # Wavelength
gaussian_size   = 1.00/sqrt(2)*mm       # 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm                 # LPForvard step size
f               = 20*cm                 # focal length of SLM lens
l               = 8                     # azimuthal order of LG beam
p               = 1
omega0_p        = .5

######## END Param ##############


####### Create LG order 
dx = (2/(N-1.))
x = y = 3*r_[-1:(1+dx):dx]
[xx, yy] = meshgrid(x,y)

lg = lg.LG_xy( l=l, p=0, xx=xx, yy=yy, omega0=2/sqrt(l) ) \
   + lg.LG_xy( l=0, p=p, xx=xx, yy=yy, omega0=.85*2/sqrt(l) )

# normalize lg; we only want the phase information
lg /= abs(lg)
#imshow (angle(lg)); show()
####### END Create LG order 



propagate = Propagate(vmin=0., vmax=100.)

F = Field( N, side_length, wavelength )
F.gaussian_aperture(gaussian_size)
#F.lens(f)
F.value[:] *= lg
#propagate(F, z=2*f, dz=.5*cm)

for z in r_[0:(2*f):(.1*f)]:
  Fcopy = F.copy()
  Fcopy.lens_forvard(f, z)
  propagate.imshow( (Fcopy.value * Fcopy.value.conj()).real )
  raw_input()
