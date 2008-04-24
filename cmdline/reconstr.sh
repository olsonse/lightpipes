#!/bin/csh -f
#	Copyright 1996, Gleb Vdovin
#	This file is a part of LightPipes distribution

#	This script imports the interferometric pattern
#	produced by interf.sh and reconstructs the
#	phase profile from imported interferogram

	
#Importing the interferogram intensity pattern:
begin 1e-2 1e-6 100 | fil_ter int subs int_in >foo

#Fourier transform and shift in the spectral domain: 
#Note the shift value (2e-3) is dependent onto the grid sampling,
#pute 1e-3 for 200x200 and 5e-4 for 400x400:

pip_fft 1 < foo | interpol 1e-2 same  2e-3 0 > foo1

#Postscript of the intensity distribution in the Fourier domain:
interpol 4e-3 < foo1 | file_ps fft_shifted.ps 100 14 |rect_ap 15e-4 1\
|file_ps fft_shifted_filtered.ps 100 14> /dev/null

#Filtering in the Fourier domain, reverse FFT and saving of the phase:
rect_ap 15e-4 1 < foo1 |pip_fft -1|file_pha pha 50 > foo

#Unwrapping the phase:
unf3 -1 < pha > pha1
#All done

