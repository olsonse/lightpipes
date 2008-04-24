@echo off

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0 distribution
rem        Another Young's interferometer        

rem model of a screen with two round holes
	begin 0.005 0.55e-6 256 > beam

rem Two circular apertures:
       circ_ap 0.00012 -0.0005 < beam > beam1
       circ_ap 0.00012 0.0005  < beam > beam2

rem Mixing two beams, rewriting initial beam:
       b_mix beam1 < beam2 > beam

rem interpolation and Postscript Picture in0r.ps
	file_ps young0.ps 128 < beam >  null

rem propagation 0.75 m and second file in1r.ps
	forvard 0.75 < beam  | file_ps young1.ps 128  > null

rem removing all beam files
	del beam 
	del beam1 
	del beam2

echo all done
