@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0 distribution
rem        Model of Young's interferometer        

echo Model of a screen with two round holes
	begin 0.005 0.55e-6 256 > beam
	rect_ap 0.0001 0.0025  -0.0005 < beam > beam1
	rect_ap 0.0001 0.0025  0.0005 0 15 < beam > beam2
rem Mixing three beams, rewriting the initial beam:
       b_mix beam1 < beam2  > beam
echo Initial field
rem show beam

rem writing  Postscript Picture 
	file_ps slit0.ps 128 < beam > null

rem propagation 0.75 m and second postscript file 
	forvard 0.75 < beam | file_ps slit1.ps 128 > beam1
echo Diffracted field
rem show beam1

rem removing all beam files
	del beam 
	del beam1 
	del beam2 

echo All done
