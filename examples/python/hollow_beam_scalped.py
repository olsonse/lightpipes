#!/usr/bin/env python
"""
LightPipes for Octave Optical Toolbox

This file is for taking a look at the optical field for doing evaporative
cooling by manufacturing a major leak at the top of the trap.  

This version looks at the result of placing the knife edge between the two
relaying lenses (that relay the SLM plane to fslm away from the MOT).
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
gaussian_size   = 0.40/sqrt(2)*mm       # 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm                 # LPForvard step size
fslm            = 15*cm                 # focal length of SLM lens
fslmR           = 30*cm                 # focal length of first relay (to move SLM plane)
fr1             = 20*cm                 # focal length of relay:2 lens:[1,4] (around MOT)
fr2             =  8*cm                 # focal length of relay:2 lens:[2,3] (around MOT)
l               = 8                     # azimuthal order of LG beam

######## END Param ##############


####### Create LG order 
dx = (2/(N-1.))
x = y = 3*r_[-1:(1+dx):dx]
[xx, yy] = meshgrid(x,y)

lg = lg.LG_xy( l=8, p=0, xx=xx, yy=yy, omega0=2/sqrt(l) )
# normalize lg; we only want the phase information
lg /= abs(lg)
#imshow (x,y,arg(lg))
####### END Create LG order 


propagate = Propagate()

if False:
  #####  Just the hollow beam relayed. 
  F = Field( N, side_length, wavelength )         # beginning field
  F.gaussian_aperture(gaussian_size)
  F.value[:] *= lg                                # SLM LG phase
  F.lens(    fslm       )                         # 'SLM' lens
  F.forvard( fslmR      )
  F.lens(    fslmR      )                         # relay:1 lens:1
  F.forvard( 2*fslmR    )
  F.lens(    fslmR      )                         # relay:1 lens:2
  F.forvard( fslmR+fslm )                         # MOT (1st pass)

  imshow( (F.value * F.value.conj()).real )
  print 'press enter to continue'
  raw_input()

  F.forvard( fr1     )
  F.lens(    fr1,    )                            # relay:2 lens:1
  F.forvard( fr1+fr2 )
  F.lens(    fr2,    )                            # relay:2 lens:2
  F.forvard( 2*fr2   )
  F.lens(    fr2,    )                            # relay:2 lens:3
  F.forvard( fr1+fr2 )
  F.lens(    fr1,    )                            # relay:2 lens:1
  F.forvard( fr1     )                            # MOT (2nd pass)
  #F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4)

  imshow( (F.value * F.value.conj()).real )



if False:
  #####  With the knife edge between the SLM (at SLM focus) relay:1,lens:1.
  F = Field(N, side_length, wavelength )          # beginning field
  F.gaussian_aperture(gaussian_size)
  F.value[:] *= lg                                # SLM LG phase
  F.lens(    fslm                        )        # 'SLM' lens
  F.forvard( fslm                        )
  F.rectangular_screen(50,100,25+100*um,0)        # knife edge
  F.forvard( fslmR-fslm                  )
  F.lens(    fslmR                       )        # relay:1 lens:1
  F.forvard( 2*fslmR                     )
  F.lens(    fslmR                       )        # relay:1 lens:2
  F.forvard( fslmR+fslm)                          # MOT (1st pass)

  ishow( (F.value * F.value.conj()).real )
  fprintf (stderr, 'press enter to continue')
  raw_input()

  F.forvard( fr1     )
  F.lens(    fr1     )                           # relay:2 lens:1
  F.forvard( fr1+fr2 )      
  F.lens(    fr2     )                           # relay:2 lens:2
  F.forvard( 2*fr2   )      
  F.lens(    fr2     )                           # relay:2 lens:3
  F.forvard( fr1+fr2 )      
  F.lens(    fr1     )                           # relay:2 lens:1
  F.forvard( fr1     )                           # MOT (2nd pass)
  #F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4)

  imshow( (F.value * F.value.conj()).real )


if True:
  #####  With the knife edge between the SLM (at SLM focus) relay:1,lens:1.
  #####  In this case, we relay the knife edge to the MOT and not the SLM plane.
  F = Field( N, side_length, wavelength )         # beginning field
  F.gaussian_aperture(gaussian_size)
  F.value[:] *= lg                                # SLM LG phase
  F.lens(    fslm           )                     # 'SLM' lens
  F.forvard( fslm           )
  F.rectangular_screen(50,100,25+100*um,0)        # knife edge
  F.forvard( fslmR          )
  F.lens(    fslmR          )                     # relay:1 lens:1
  F.forvard( 2*fslmR        )
  F.lens(    fslmR          )                     # relay:1 lens:2
  F.forvard( fslmR          )                     # MOT (1st pass)

  imshow( (F.value * F.value.conj()).real )
  print 'press enter to continue'
  raw_input()

  F.forvard( fr1     )
  F.lens(    fr1     )                            # relay:2 lens:1
  F.forvard( fr1+fr2 )
  F.lens(    fr2     )                            # relay:2 lens:2
  F.forvard( 2*fr2   )
  F.lens(    fr2     )                            # relay:2 lens:3
  F.forvard( fr1+fr2 )
  F.lens(    fr1     )                            # relay:2 lens:1
  F.forvard( fr1     )                            # MOT (2nd pass)
  #F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4)

  imshow( (F.value * F.value.conj()).real )

print 'press enter to exit'
raw_input()
