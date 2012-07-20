#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox

creates a Laguerre-Gaussian beam (actually a non pure LG mode by projecting
a Gaussian onto LG(l=1)) and propagates this through a lens, grating, piece
of glass 1mm thick, to a distance a little further away than the lens focal
length.
"""

from common import Propagate
from lightpipes import Field, lg, grating
from pylab import *

m=1
nm=1e-9*m
um=1e-6*m
mm=1e-3*m
cm=1e-2*m
W=1
mW=1e-3*W

N               = 512                   # Number of pixels
grid_dx         = 15.*um                # grid cell size
side_length     = N* grid_dx            # physical size of SLM
wavelength      = 780*nm                # Wavelength
gaussian_size   = 1.00/sqrt(2)*mm       # 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm                 # LPForvard step size
f               = 20*cm                 # focal length of SLM lens

lf              = 0.64                  # fraction of f to travel before going through glass.
r               = 0.05                  # reflection at each surface of glass.
a0              = 45*pi / 180           # incident angle of primary beam
n               = 1.55                  # index of refraction of glass
a               = sin(a0)/n             # angle of beams in glass (Snell's law)
k               = n*2*pi/wavelength     # k in glass
d               = 1.*mm / cos(a)        # 1-pass distance travelled through glass of secondary beam
dphi            = 2*k*d                 # relative phase delay on the secondary beam
dx              = 2 * d*sin(a)          # displacement of secondary beam parallel to glass.
dx0             = dx * cos(a0)          # displacement of secondary beam transverse to output beam.
l               = 8                     # azimuthal order of LG beam

P               = 250*mW                # beam power

######## END Param ##############


####### Create LG order 
dx = (2/(N-1.))
x = y = (side_length/2)*r_[-1:(1+dx):dx]
[xx, yy] = meshgrid(x,y)

lg = lg.LG_xy( l=l, p=0,
               xx=xx*6/side_length,
               yy=yy*6/side_length,
               omega0=2/sqrt(l) )
# normalize lg; we only want the phase information
lg /= abs(lg)
#imshow (x,y,angle(lg))
####### END Create LG order 



propagate = Propagate()

F = Field( N, side_length, wavelength )
F.gaussian_aperture(gaussian_size)

Pf   = P/sum(F.value * F.value.conj())  # power per grid cell [ W ]
Ef   = sqrt(Pf/grid_dx**2)              # sqrt(I) [ sqrt(W/m^2) ]
Ef  /= sqrt(mW/cm**2)                   # put Ef in units [ sqrt(mW/cm^2) ]
F.value[:] *= Ef                        # put F in units of sqrt(I)

print 'writing SLM phase...'
F.lens(f)             # apply lens (physical lens)
Fb = F.copy()         # store the undeflected beam
Fb.value[:] *= sqrt(0.2) # undeflected light
F.value[:] *= lg      # apply lg phase
grating(F,32,0)       # apply a grating phase.
F.value[:] *= sqrt(0.8) # deflection efficiency of SLM
F.value[:] += Fb.value # mix the beams back together.
print 'finished writing SLM phase.'

F.forvard(lf*f)

# traverse a piece of glass and interfere primary beam with 1st internally
# reflected+transmitted beam.
Fr = F.copy()
Fr.interpolate(side_length,N,dx0)
Fr.value[:] *= exp(1j * dphi)     # add phase delay to secondary beam
Fr.value[:] *= (r*r*(1-r))        # reflection, reflection, transmission
F.value[:]  *= (1-r)              # transmission
F.value[:]  += Fr.value           # mix the beams back together.


xlabel('x (mm)')
ylabel('y (mm)')
propagate(F, z=(1-lf)*f + 2*cm, dz=.5*cm)
cb = colorbar()
cb.set_label('I [mW/cm^2]')
print 'press enter to exit'
raw_input()


#  results:  
# Pure LG mode (+grating+lens only), the darkness is down
# to:  ~0.00001 mW/cm^2
#
# Gaussian beam projected onto LG(l=1) mode:
# with neither the non-diffracted fraction of light nor the secondary beam
# from the glass transmission, the darkness at the center of the trap is down
# to   ~0.00185 mW/cm^2

# Adding the secondary beam tranmission from the glass, the darkness rises to:
#      ~0.0300 mW/cm^2

# Adding also the portion of non-diffracted light (assuming the physical
# lens), the darkness is:
#      ~0.022 mW/cm^2 (must have destructively interfered at the center

