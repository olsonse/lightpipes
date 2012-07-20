#!/usr/bin/env python
"""
LightPipes for Octave Optical Toolbox
"""

from common import Propagate
from lightpipes import Field, lg
from pylab import *

m=1
nm=1e-9*m
um=1e-6*m
mm=1e-3*m
cm=1e-2*m

N               = 512           # Number of pixels
grid_dx         = 15*um         # pixel size on slm
side_length     = N * grid_dx   # physical size of SLM
wavelength      = 780*nm        # Wavelength
f_torus         = 40*cm         # lens applied directly on camera
R0              = .25*N*grid_dx # radius of center of toroidal lens
R1              = .19*N*grid_dx # central radius of phase ring.
w1              = .23*N*grid_dx # width of phase ring.
#w               = 0.12*40*mm/sqrt(2) # 40mm colimating lens
w               = 1.71*mm/sqrt(2) # smaller beam

def tlens(x,y,R0,f,wavelength):
  """
  This function calculates the phase to use for a
  toroidal lens + quadratic lens.
  """
  return 2*pi + \
    (-pi*(sqrt(x**2+y**2) - R0)**2 / (f * wavelength) ) % (2*pi)


def phase_ring(x, y, R, w):
  """
  This function calculates the phase ring that is to be added to the phase of
  the toroidal lens at the SLM.
  """
  phi = zeros(x.shape)
  r2 = x**2 + y**2
  Ri_2 = (R - 0.5*w)**2
  Rf_2 = (R + 0.5*w)**2
  phi[ nonzero( r2 < Rf_2 ) ] = pi
  phi[ nonzero( r2 < Ri_2 ) ] = 0
  return phi


x = y = grid_dx * r_[-255.5:256:1.]
[xx, yy] = meshgrid(x,y)
torus_phase = tlens(xx,yy,R0,f_torus,wavelength) + phase_ring(xx,yy,R1,w1)

propagate = Propagate()

#print 'press enter to continue'
#propagate.imshow( torus_phase )
#raw_input()

# propagate the field from the SLM and take a look
F = Field( N, side_length, wavelength)     # initial field
F.gaussian_aperture(w)                     # gaussian input
F.circular_aperture(side_length/2.)        # Iris in front of SLM
F.value[:] *= exp(1j * torus_phase )       # apply SLM phase
F.forvard(f_torus)
propagate(F, z=f_torus, dz=f_torus/25.)

print 'press enter to show final phase...'
raw_input()
propagate.imshow( angle(F.value) )

print 'press enter to show final intensity...'
raw_input()
propagate.imshow( (F.value * F.value.conj()).real )

print 'press enter to exit'
raw_input()
