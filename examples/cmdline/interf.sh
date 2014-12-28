#!/bin/sh
#	Copyright 1996, Gleb Vdovin,
#	This file is a part of LightPipes distribution

# 	Simple interferogram is formed here, to
#	be imported and reconstructed by script
#	reconstr.sh
#	Note: the sampling is important for
#	further processing, if sampling is changed,
#	the value of frequency shift to be changed in
#	reconstr.sh.


#	Aberrated beam - you may try different aberrations here:
begin 1e-2 1e-6  100 | zernike 7 1 0.72e-2 10 > foo

#	Saving and unfolding the phase:
file_pha pha 50 < foo > /dev/null
unf3 1 < pha > pha1  

#	Addition of tilt, to obtain the vertical fringes:
tilt 2e-3 0 < foo > foo1

#	Ideal beam, mixing -> generation of interferogram:
begin 1e-2 1e-6  100 | b_mix foo1 | file_int int_in  > foo2

#	Plot the interferogram:
file_ps interf.ps 100 1 < foo2 > /dev/null

#	All done
