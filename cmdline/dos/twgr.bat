@echo off
break on
rem     Twyman-Green interferometer
rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0  distribution
rem        


begin 0.015 5e-7 150 | circ_ap 0.005 | Fresnel 0.5 | b_split foo 0.3 > foo1
rem  Propagation to mirror1, reflection,
rem  propagation back and beamsplitter again
Fresnel 1 <foo | absorber  0.7 >foo2
rem  Propagation to mirror2, reflection,
rem  propagation back and beamsplitter again
Fresnel 0.4 < foo1 | Zernike 3 1 0.005 25 | Fresnel 0.4 | absorber  0.3 > foo3
rem  mixing the two beams and propagation to the screen
b_mix foo2 < foo3 | Fresnel 1 | file_ps tg.ps > foo1

rem show foo1

rem  deleting the temp files
del foo 
del foo1 
del foo2 
del foo3
