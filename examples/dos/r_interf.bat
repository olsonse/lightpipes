@echo off

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0 distribution
rem        
rem        rotational interferometers:


rem        coma 
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 3 1 0.01 10 | b_split foo1 > foo
rem  rotating and mixing
interpol 0.04 256 0. 0. 180 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps rt1.ps 128 < foo2 > null


rem        astigmatism
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 2 2 0.01 10 | b_split foo1 > foo
rem  rotating and mixing
interpol 0.04 256 0. 0. 90 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps rt2.ps 256 < foo2 > null



rem        high order
begin 0.04 5e-7 | circ_ap 0.01 | Zernike 7 3 0.01 10 | b_split foo1 > foo
rem  rotating and mixing
interpol 0.04 256 0. 0. 90 < foo1 | b_mix foo > foo2
rem  postscript output
file_ps rt3.ps 128 < foo2 >null
del foo 
del foo1 
del foo2
