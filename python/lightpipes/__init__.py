import _lightpipes as lib
from _lightpipes import *
import randomize
from grating import grating
import lg


# proxy _lightpipes.Field constructor

class Field(lib.Field):
  def __init__(self, number, side_length, *a, **kw):
    if type(number) is not lib.SizePair:
      try:
        number = lib.SizePair(*number)
      except:
        number = lib.SizePair(number)
    if type(side_length) is not lib.DoublePair:
      try:
        side_length = lib.DoublePair(*side_length)
      except:
        side_length = lib.DoublePair(side_length)
    super(Field,self).__init__( number, side_length, *a, **kw )

  __init__.__doc__ = lib.Field.__init__.__doc__
