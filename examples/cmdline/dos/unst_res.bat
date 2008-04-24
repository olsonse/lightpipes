@echo off
break on


rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0  distribution
rem        

rem        LightPipes script: a model of an unstable resonator
rem        with N=10, M=2, 
rem        See D.B. Rench, Applied  Optics 13, 2546...2561 (1974)
rem        All dimensions are taken from the table on p. 2552
rem        in this reference

rem        theoretically the initial field distribution may be random 
rem        but we'll use a good plane wave because it will take too long
rem        to converge starting with a random distribution
	
begin 30. 3e-4 512 > field1

rem        40 iterations inside the resonator
rem **** copy this 40 times *********************
rem        reflection from the convex mirror and propagation 
rect_ap 5.48 < field1 | lens -1e4 | forvard 1e4 > field2
rem        reflection from the concave mirror and propagation
rect_ap 10.96 < field2 | lens 2e4  | forvard 1e4 | normal y > field1
Strehl y < field1 > null
file_int in < field1 | cros_out cross > null
rem ******* end of loop block ********************

rem ********* repeat here 20 times ***************

rem  output of the resonator is screened by the output mirror:
rem  the output 
rect_screen 5.48 < field1 | file_int in | file_pha pha | file_ps res_out.ps | Strehl y > null
rem  removing the temp files
del field1 
del field2
del null







