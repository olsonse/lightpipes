@echo off

rem       (C) Gleb Vdovin 1993-1995
rem       This file is a part of LightPipes.1.0 distribution
rem       
rem       Far field distributions for aberrated beams:
rem       Here the trick with coordinate transfer is used:
rem       first we define the field in a "small" region
rem       then we apply a weak positive lens
rem       then we apply coordinate transfer equivalent to a small
rem       negative lens (resulting in a plane wave in divergent coordinates
rem       and then we propagate it. As the coordinates are 
rem       divergent, the output field is wide enough to
rem       to match the diffracted pattern
rem       Without this trick the input grid would be too coarse

echo Ideal beam:
begin 0.1 1e-6 150 | circ_ap 0.01 |Zernike 3 1 0.01 0 | lens 200 > foo
lens_fresn -200 1000 < foo  | file_ps f_fi1.ps 150 2.5 > foo1

echo Coma with amplitude 6.28 rad or ~ 1 lambda:
begin 0.1 1e-6 150 | circ_ap 0.01 |Zernike 3 1 0.01 3.14 | lens 200 > foo
lens_fresn -200 1000 < foo | file_ps f_fi2.ps 150 2.5 > foo1

echo Astigmatism with amplitude 6.28 rad or ~ 1 lambda:
begin 0.1 1e-6 150 | circ_ap 0.01 |Zernike 2 2 0.01 3.14 |lens 200 > foo
lens_fresn -200 1000 < foo  | file_ps f_fi3.ps 150 2.5 > foo1


echo Spherical aberration with amplitude 6.28 rad or ~ 1 lambda:
begin 0.1 1e-6 150 | circ_ap 0.01 | Zernike 4 0 0.01 3.14 | lens 200 > foo
lens_fresn -200 1000 < foo  | file_ps f_fi4.ps 150 2.5 > foo1

echo Compare this with  direct calculation
echo without the coordinate transform:
begin 0.6 1e-6 150 | circ_ap 0.01 | Zernike 4 0 0.01 3.14 > foo
Fresnel  1000 < foo  | file_ps f_fi4_b.ps 150 2.5 > foo1
del foo
del foo1


