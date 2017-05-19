"""
LightPipes for Octave Optical Toolbox
"""

import numpy as np

#################### BEGIN GRATING FUNCTIONS ###############################
def perfect(F,px,py):
  """
  applies a perfect grating phase with px,py periods in x,y directions
  respectively.  This is indicated as perfect since it will directly apply a
  periodic ramp up to 2*pi phase shift, regardless of wavelength.
  """
  Gx = make_perfect_grating_phase(F.info.number,px)
  Gy = make_perfect_grating_phase(F.info.number,py)

  for k in xrange(F.info.number):
    F.value[:,k] *= Gx[k]
    F.value[k,:] *= Gy[k]
  return F



def make_perfect_grating_phase(N,p):
  """
  returns the field of a grating function with 'p' periods over the N pixels.
  """
  return np.exp(1j * np.r_[0:N] * ( 2*np.pi * p / N ) )



class Base(object):
  def G(self, number, side_length):
    """
    return the phase array
    """
    return np.exp(1j * self(number, side_length))

  def apply(self, F, axis=0):
    if axis == 0:
      F.value[:] = (
        F.value.transpose() * self.G(F.info.number.first, F.info.side_length.first)
      ).transpose()

    elif axis == 1:
      F.value[:] *= self.G(F.info.number.second, F.info.side_length.second)

    else:
      raise RuntimeError('Can only apply grating in either x or y direction')



class Blazed(Base):
  def __init__(self, period, wavelength, design_wavelength=None, turning_point=1):
    """
    turning_point : The grating does not jump immediately to the
      beginning of the next saw-tooth, but rather turns back down as soon as the
      total phase (for the design wavelength) has reached
      (turning_point * (2*pi)).
    """
    super(Blazed,self).__init__()
    if design_wavelength is None:
      design_wavelength = wavelength

    self.period = period
    self.wavelength = wavelength
    self.design_wavelength = design_wavelength

    self.k = 2*np.pi / self.wavelength
    # blaze angle, assuming that the grating is at "littrow configuration" for the
    # design_wavelength
    self.blaze = np.arcsin(design_wavelength / (2*period))

    # the distance it takes for the phase to change by 2*pi when operated at the
    # design_wavelength
    self.period_prime = turning_point \
                      * self.design_wavelength / (2*np.tan(self.blaze))


  def __call__(self, number, side_length):
    x = (np.r_[0:number] * side_length / number) % self.period
    retval = (2*self.k * np.tan(self.blaze) * x)

    # now take care of the grating sloping back down to the bottom again
    i = np.argwhere(x >= self.period_prime)
    retval[i] = self.k * self.design_wavelength \
              * (self.period - x[i]) \
              / (self.period - self.period_prime)
    return retval





class Binary(Base):
  """
  Creates a binary grating that shifts phase by depth.
  """
  def __init__(self, period, depth=np.pi, duty=0.5):
    super(Binary,self).__init__()
    self.period = period
    self.depth  = depth
    self.duty   = duty

  def __call__(self, number, side_length):
    dx = side_length / (number - 1)
    return self.depth * (((dx * np.r_[0:number]) % self.period) < (self.duty*self.period))
###################### END GRATING FUNCTIONS ###############################
