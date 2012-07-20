# vim: ts=2:sw=2:tw=80:nowrap
import sys
sys.path.insert(0,'../../python') # prepend this to get this version

import numpy as np
import pylab

class Propagate:
  def __init__(self):
    self.im = None
    pylab.ion()

  def __call__(self, F, z, dz):
    if self.im:
      self.im.set_array( (F.value * F.value.conj()).real )
    else:
      self.im = pylab.imshow( (F.value * F.value.conj()).real )

    for zi in np.arange(0, z+dz, dz):
      F.forvard(dz)
      self.im.set_array( (F.value * F.value.conj()).real )
      pylab.draw()
