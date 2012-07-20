#!/usr/bin/env python
"""
LightPipes for Python Optical Toolbox
"""

import common
from lightpipes import Field
from pylab import *


phi = 178 * (pi/180.)


F = Field(256, 512*15e-6, 780e-9)
F.gaussian_aperture(.75e-3) 
F.axicon(phi, 1.5,   0e-6, 0)
F.forvard(25e-2)

imshow( abs(F.value)**2 )
show()

F.axicon(phi, 1.5,   0e-6, 0)
F.forvard(100e-2)

imshow( abs(F.value)**2 )
show()


#F.lens(-20e-2)
#F.forvard(3e-2)
#    gsplot F7.F.*conj(F7.F);
