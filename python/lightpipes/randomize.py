from numpy import exp
from numpy.random import rand

def intensity(F):
  """
  Substitutes a random intensity in the field.
  This random intensity is normalized to 1.  Use Fout *= 1e-3 to
  normalize to 1 milliwatt for example.
  """
  F.value[:] = F.value / abs(F.value) * rand(*F.value.shape)


def phase(mx,F):
  """
  Substitutes a random phase in the field.
    mx : maximum value of the random phase
  """
  F.value[:] = abs(F.value) * exp(1j * mx * rand(*F.value.shape) )
