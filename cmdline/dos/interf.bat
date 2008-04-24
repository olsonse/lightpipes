@echo off
rem     This script generates interferogram which
rem     is later reconstructed by reconstr.bat

rem        Aberrated beam:
begin 1e-2 1e-6 200 | Zernike 7 1 0.72e-2 10 > foo

rem        Saving and unfolding the phase:
file_pha pha 50 < foo > null
unf3 1 < pha > pha1  

rem        Addition of tilt, to obtain the vertical fringes:
tilt 2e-3 0 < foo > foo1

rem        Ideal beam, mixing -> generation of interferogram:
begin 1e-2 1e-6 200 | b_mix foo1 | file_int int_in 256 > foo2

rem        Plot the interferogram:
file_ps interf.ps 200 1 < foo2 > null
del null
rem        All done
