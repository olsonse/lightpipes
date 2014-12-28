@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes beta distribution
rem        

echo Shearing intererometer, defocus
echo starting, lens, Z1 

begin 0.04 5e-7 128 | circ_ap 0.01 | lens -20 | forvard 0.5 | b_split foo1 > foo
rem  shifting beam 3 mm without real interpolation, mixing it
rem  with centered beam
interpol 0.04 128 0.003 0.001 < foo1 | b_mix foo | forvard 0.5 > foo2
rem  postscript output
file_ps shear1.ps 128 < foo2 | cros_out out > null
rem show foo2


echo Shearing intererometer, astigmatism
begin 0.04 5e-7 128 | circ_ap 0.01 | Zernike 2 2 0.01 20 | forvard 0.5 | b_split foo1 > foo
rem  shifting beam 3 mm without real interpolation, mixing it
rem  with centered beam
interpol 0.04 128 0.003 0.0 < foo1 | b_mix foo | forvard 0.5 > foo2
rem  postscript output
file_ps shear5.ps 128 < foo2 > null
rem show foo2


echo Shearing intererometer, coma
begin 0.04 5e-7 128 | circ_ap 0.01 | Zernike 3 1 0.01 10 | forvard 0.5 | b_split foo1 > foo
rem  shifting beam 3 mm without real interpolation, mixing it
rem  with centered beam
interpol 0.04 128 0.003 0.00 < foo1 | b_mix foo | forvard 0.5 > foo2
rem  postscript output
file_ps shear4.ps 128 < foo2 > null
rem show foo2


echo This may be a spherical aberration ??
begin 0.04 5e-7 128 | circ_ap 0.01 | Zernike 4 0 0.01 10 | forvard 0.5 | b_split foo1 > foo
rem  shifting beam 3 mm without real interpolation, mixing it
rem  with centered beam
interpol 0.04 128 0.003 0.00 < foo1 | b_mix foo | forvard 0.5 > foo2
rem  postscript output
file_ps shear3.ps 128 < foo2 > null
rem show foo2
del  foo 
del foo1 
del foo2
del null


