"""
LightPipes for Octave Optical Toolbox
"""

import numpy as np

#################### BEGIN GRATING FUNCTIONS ###############################
def grating(F,px,py):
  """
  applies a grating with px,py periods in x,y directions respectively.
  """
  Gx = make_grating(F.info.number,px)
  Gy = make_grating(F.info.number,py)

  for k in xrange(F.info.number):
    F.value[:,k] *= Gx[k]
    F.value[k,:] *= Gy[k]
  return F



def make_grating(N,p):
  """
  returns the field of a grating function with 'p' periods over the N pixels.
  """
  return np.exp(1j * np.r_[0:N] * ( 2*np.pi * p / N ) )



###################### END GRATING FUNCTIONS ###############################
