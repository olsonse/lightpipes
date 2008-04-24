@echo off

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0 distribution
rem        model of a radial shearing interferometer        


rem defocus
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 2 0 0.01 10 | b_split foo1 > foo
rem  Radial magnification
interpol 0.04 256 0. 0. 0 1.3 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps radi1.ps 128 < foo2  > null
rem  This will work if you have ghostview installed, otherwise comment it
rem show foo2


rem spherical aberration
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 4 0 0.01 10 | b_split foo1 > foo
rem  Radial magnification
interpol 0.04 256 0. 0. 0 1.3 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps radi2.ps 128 < foo2 > null
rem  This will work if you have ghostview installed, otherwise comment it
rem show foo2


rem High order  aberration without axial symmetry
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 10 4 0.01 10 | b_split foo1 > foo
rem  Radial magnification
interpol 0.04 256 0. 0. 0 1.3 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps radi3.ps 128 < foo2 > null
rem  This will work if you have ghostview installed, otherwise comment it
rem show foo2
del foo 
del foo1 
del foo2
