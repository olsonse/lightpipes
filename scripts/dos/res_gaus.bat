@echo off

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a LightPipes 1.0 distribution
rem        
rem -------------------------------------------------
rem        Model of a solid-state laser with a thermal lens
rem        Base=1m
rem        ouput mirror has Gaussian reflection
rem        with Rmax =60%
rem        plane parallel mirrors
rem        Thermal lens with F=5m
rem        active rod has a diameter of 6mm  
	
set num_point = 120
set num_iter= 5

rem  here we form a pseudo-random phase:

begin 12  1e-3 $num_point | random 6.28 > field1 

echo  $num_iter  iterations inside the resonator
echo  (in a demo I put only 1 iteration)

rem ----------- copy this block $num_iter times ------------------------
gauss 2.5  0 0 0.6 < field1 | Fresnel 500 | circ_ap 3 | lens 5000 > field2 
Fresnel 500 < field2 | lens -2500 | Fresnel 500 | lens 5000 > foo
circ_ap 3 < foo |l_amplif 1e-2 100 0.5 | Fresnel 500 | file_ps ps > field1 
gauss_screen 2.5 0 0 0.6 < field1 | Strehl y > null 
rem --------------------end of block ------------------------------------

echo End
rem  output of the resonator is screened by the output mirror:
echo  The output 
gauss_screen 2.5 0 0 0.6 < field1 | file_int in | file_pha pha | file_ps ps > null

