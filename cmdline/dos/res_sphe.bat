@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a LightPipes 1.0 distribution
rem        Nodel of a laser with unstable resonator 
rem        and active medium with gain 1e-4 mm^{-1}
rem        length of 1e4 mm and I_sat of 1W/(mm^2)
rem 
rem        all dimensions are in mm here

	
begin 7. 3e-4 100 > field1

rem ------   copy this block $num_iter times (40 is good)    -----
rem        40 iterations inside the resonator
rem        reflection from the convex mirror and propagation 
rect_ap 5.48 < field1 | l_amplif 1e-4 1e4 1 | lens_fresn -1e4  1e4 > field2
rem        reflection from the concave mirror and propagation
rect_ap 10.96 < field2 | l_amplif 1e-4 1e4 1 | lens_fresn 2e4  1e4 | Strehl y > field1
copy field1 field_out
interpol 7. < field1 > field2 
move field2 field1
rem -------------------end of the loop block----------------
rem ------------don't forget to repeat 40 times ;-) -----------------

rem  output of the resonator is screened by the output mirror:
rem  the output 
echo Output beam:
convert y < field_out | rect_screen 5.48 | file_ps res_out.ps | Strehl y > field1
rem before screening with the output mirror
rem show field_out
rem after screning with the output mirror
rem show field1

rem  removing the temp files
del field1 
del field_out
