@echo off
rem Reconstruction of the interferogram (created by interf.bat):
rem Importing the interferogram intensity pattern:
begin 1e-2 1e-6 200 | fil_ter int subs int_in > foo

rem Fourier transform and shift in the spectral domain: 
pip_fft 1 < foo | interpol 1e-2 200  10.e-4 0 > foo1

rem Postscript of the intensity distribution in the Fourier domain:
interpol 4e-3 < foo1 | file_ps fft_shi.ps 200 14 | rect_ap 15e-4 1 | file_ps fft_sf.ps 200 14 > null

rem Filtering in the Fourier domain, reverse FFT and saving the phase:
rect_ap 15e-4 1 < foo1 | pip_fft -1 | file_pha pha 50 > foo

rem Unwrapping the phase:
unf3 -1 < pha > pha1

rem All done

