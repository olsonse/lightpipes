@echo off
break on

rem       (C) Gleb Vdovin 1993-1995
rem       This file is a part of LightPipes.1.0
rem       LightPipes.1.0 script: demonstration of amplitude and 
rem       phase filtering

rem       Here a shifted Gaussian distribution 
rem       aberrated with coma is formed   

begin 0.1 1e-6 256 | gauss 0.02 0.02 0.01 | zernike 3 1 0.02 1 > field1

rem       Intensity and phase are saved in files in and pha
rem       Note the command line argument 256 (!) -
rem       the default is 128

file_int in 256 < field1 |file_pha pha 256 > null

rem       Here we form a new beam with a uniform 
rem       phase and intensity:

begin 0.1 1e-6 256 > field2

rem       Here we filter it through the saved intensity and phase filters:
rem       The field in file field3 has Gaussian intensity 
rem       and uniform phase

fil_ter int subs in  < field2 > field3 

rem       The field in file field4 has uniform intensity
rem       and aberrated phase

fil_ter pha subst pha < field2 > field4

rem       The field in a file field5 has the same distribution
rem       of the complex amplitude as the initial beam field1

fil_ter pha subs pha <field2 | fil_ter int subs in > field5

rem       End, all the files have to be erased:

del field1 
del field2 
del field3 
del field4 
del field5
del null





