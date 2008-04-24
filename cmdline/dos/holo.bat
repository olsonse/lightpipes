@echo off
break on

rem        (C) Gleb Vdovin 1993-1995
rem        This file is a part of LightPipes.1.0  distribution
rem        
rem        synthesis of a Fresnel lens:

rem spherical wavefront -> beam

begin 0.015 5e-7 256 | circ_ap 0.005 | lens 2.5 > foo1

rem a beam with a plane wavefront

begin 0.015 5e-7 256 | circ_ap 0.005  > foo2
b_mix foo1 < foo2 | file_ps holo.ps > foo3

rem show the pattern:
rem show foo3

rem  deleting the temp files
del foo1 
del foo2 
del foo3
