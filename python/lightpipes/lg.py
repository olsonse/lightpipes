# vim: ts=2:sw=2:tw=80:nowrap
"""
LightPipes for Python Optical Toolbox

Provision for Laguerre-Gaussian beams.
"""

from numpy import *

def fact(n):
  """
  returns n! where n is an integer
  """
  n = floor(n)
  if (n > 0):
    #retval = n * fact (n-1);

    retval = 1
    while n >= 1:
      retval *= n;
      n -= 1
    return retval

  return 1



def LG(l,p,rho,phi,omega0):
  """
  returns the field of a laguerre gaussian beam of order (l,p).  Note
  that the result is unitless (i.e. scalable in position space).
  where:
          l   = azimuthal order
          p   = radial order
          rho = radial position of evaluation of LG mode
          phi = azimuthal position of evaluation of LG mode
          omega0 = 
  Note that the electric field amplitude is given by
      (sqrt(P_{0})/omega0) * LG(l,p,rho,phi,omega0)
  where P_{0} is the beam power.
  """
  return sqrt(2 * fact(p) / (pi * fact(p+abs(l))))      \
         * (sqrt(2) * rho / omega0)**(abs(l))           \
         * LaguerreL( abs(l), p, 2*rho**2 / omega0**2 ) \
         * exp( -rho**2/omega0**2 + 1j*l*phi)


def PHI(x,y):
  """
  return the positive (CCW) angle from the positive x axis.  This
  function correctly determines the angle independent of quadrant.
  where:
          x = x position
          y = y position
  """
  ANGLE = acos(x/sqrt(x**2 + y**2))
  k = find(y < 0)
  ANGLE[k] = 2*pi - ANGLE(k)
  return ANGLE


def LG_xy(l,p,xx,yy,omega0):
  """
  l   = azimuthal order
  p   = radial order
  xx  = x position of evaluation of LG mode
  yy  = y position of evaluation of LG mode
  omega0 = 
  Note that this function is defined in terms of LG(...) and PHI(x,y)
  """
  return LG(l,p,sqrt(xx**2+yy**2), angle(xx + yy*1j), omega0)


def LaguerreL(l,p,rho):
  """
  returns the general Laguerre function L(l,p,rho)
  where:
          l = azimuthal order
          p = radial order
          rho = radial position
  """
  genL = 0
  for m in xrange(0,p+1):
    genL += (-1)**m / (fact(p-m) * fact(abs(l)+m) * fact(m)) * rho**m

  genL *= fact(abs(l)+p)
  return genL
