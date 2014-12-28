@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes distribution
rem        Here we form two intensity distributions        
rem        and then reconstruct the phase in the near and far field
rem        using these distributions

	
set num_point=200
set num_iter=5
set dist = 1

echo Initial distribution is formed here:
echo a model of a screen with two slits
       
	begin 0.005 0.25e-6 $num_point > beam
       rect_ap 0.0001 0.0025  -0.0005 < beam > beam1
       rect_ap 0.0001 0.0025  0.0005 0 15 < beam > beam2
echo Mixing the beams, rewriting the initial beam:
       b_mix beam1 < beam2  > field_near

file_int int_near $num_point < field_near > null



echo The result of diffraction is here:

forvard 0.2 < field_near > field_far

file_int int_far $num_point < field_far > null


begin 0.005 0.25e-6 $num_point | fil_ter int subs int_near > ref_file

echo Here a random field is formed to start with 

random 6.28 < field_near | forvard 0.2 > f_far1

echo 5  iterations to reconstruct the phase


	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       
	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       
	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       
	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       
	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       
	fil_ter int subs int_far < f_far1 | forvard -0.2 > f_near1
	fil_ter int subs int_near < f_near1 | forvard 0.2 > f_far1       


rem        use "show" to see how f_near1 and f_far1 look like
	
rem show f_far1
rem show f_near1
del null










