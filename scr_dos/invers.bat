@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0 distribution
rem        
rem 
rem        this script demonstrates the possibility of
rem        phase reconstruction from a two near-field intensity
rem        distributions

echo We start with a plane phase:

begin  0.011 6.33e-7 128 > f_far1


rem        num_iter  iterations 


rem ------iterated block, repeat 100 times ;-) ------------------
echo       Here we substitute the known second intensity and interpolate
echo       into a bigger grid - to make "forvard" happy
	fil_ter int subs i_2 < f_far1 | interpol 0.02 256 > foo

echo       Propagation "back" to the plane of the first image
	forvard -2 < foo |interpol 1.1e-2 128 > f_near1

echo       substitution of the first intensity distribution
echo       and propagation forward:
	fil_ter int subs  i_1 < f_near1 | Fresnel 2 > f_far1
rem ------------end of iterated block ---------------------------


rem ------------copy it to here 99 times ------------------------

echo After all the phase of f_near1 is the reconstructed one.




 







