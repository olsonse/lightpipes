#!/usr/bin/env python
"""
LightPipes for Octave Optical Toolbox

This file demonstrates diffraction from a binary grating where the optical field
is not square, thus saving computational resources while still resolving the
grating.
"""

from physical.unit import *
from common import Propagate1D
from lightpipes import Field, grating
from pylab import *


def main(return_final=False):
  N               = r_[2**21,2]           # Number of pixels (different resolution for demo)
  DX              = 60*mm                 # physical size in the important direction
  wavelength      = 780*nm                # Wavelength
  gaussian_size   = 1.0/sqrt(2)*mm       # 1/e Intensity width of beam; (note, it's not the 1/e^2)
  L_propagate     = 25*mm                 # distance from initial lens to 2D-MOT

  # grating parameters
  grating_period  = mm/900
  grating_depth  = np.pi                  # for the binary grating


  ######## END Param ##############

  side_length     = DX * r_[1, N[1]/float(N[0])] # physical size


  ####### Create grating
  G = grating.phase.Binary(grating_period/m, grating_depth)
  ####### END Create grating


  propagate = Propagate1D(N, side_length, window_size=511)

  #####  With the knife edge between the SLM (at SLM focus) relay:1,lens:1.
  #####  In this case, we relay the knife edge to the MOT and not the SLM plane.
  F = Field( N, side_length/m, wavelength/m )   # beginning field
  F.gaussian_aperture(gaussian_size/m, 1e300)

  G.apply(F)

  if return_final:
    F.forvard(L_propagate/m)
    plot( (F.value * F.value.conj()).real[:,N[1]/2] )
    return F
  else:
    propagate(F, L_propagate/m, L_propagate/m/75.)


  print 'press enter to exit'
  raw_input()

if __name__ == '__main__':
  main()
