#!/usr/bin/env python
"""
LightPipes for Octave Optical Toolbox
Simulates a unstable resonator
"""

import common
from lightpipes import Field, randomize
from pylab import *
import time

m=1
nm=1e-9*m
mm=1e-3*m
cm=1e-2*m

wavelength  = 308*nm
size        = 14*mm
N           = 100
w           = 5.48*mm
f1          = -10*m
f2          = 20*m
L           = 10*m
Isat        = 1.0
alpha       = 1e-4
Lgain       = 1e4

F = Field(N,size,wavelength)
randomize.intensity(F)
randomize.phase(1,F)

SR = dict()
for l in xrange(10):
  F.rectangular_aperture(w,w)
  F.l_amplify(alpha,Lgain,Isat)
  F.lens_fresnel(f1,L)
  F.l_amplify(alpha,Lgain,Isat)
  F.lens_fresnel(f2,L)
  SR[l]=F.get_strehl()
  F.interpolate(size,N)
  print 'Round trip {l} Strehl ratio= {s}'.format(l=l,s=SR[l])

  F2 = F.copy()
  F2.rectangular_screen(w,w)

  imshow( abs(F2.value)**2 ) # plot intensity
  show()
  time.sleep(.3)

F2.spherical_to_normal_coords()

print 'SR: ', SR
SR = SR.items()
SR.sort( key = lambda v : v[0] )
SR = array( SR )
figure(2)
plot(SR[:,0], SR[:,1])
xlabel('Number of Roundtrips')
ylabel('Strehl ratio')

figure(3)
imshow( abs(F2.value)**2 ) # plot intensity
xlabel(''); ylabel('')
title('Intensity distribution just behind the outcoupler')

# Far-field calculation:
z=1*m
f=40*m
ff = z*f / float(f-z)
F2.lens(f)
F2.lens_fresnel(ff,z)
F2.spherical_to_normal_coords()

figure(4);
imshow( abs(F2.value)**2 ) # plot intensity
xlabel(''); ylabel('');
title('Intensity distribution in the far field')
show()
