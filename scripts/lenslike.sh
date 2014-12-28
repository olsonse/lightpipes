#!/bin/csh -f

#	We form here the distribution of the refractive index
#	in the file refr.ind;
#       med is a Fortran program to do that:
echo "You need a fortran compiler to run this script"

f77 -o med  med.f -lm 
echo 100 1.5 400 1e-3  | med > refr.ind

#	The array of refractive indexes (1.5 ) is used here
#	to form an array of absorbtion coefficients.
#	No optics here, it's just a reuse of the table of numbers!
#	 

scale -1. < refr.ind > abs.co

# 	propagation of a  gauss beam in a quadratic medium: 
#	the steps back are implemented with the 
#	second operator to show the reversibility of the
#	operator:


begin 1e-3 1e-6 100 |  gauss 1.13e-4 | Strehl y | tilt 1e-3 0|\
steps 4e-3 250  refr.ind  abs.co s.out 5 | Strehl y |\
steps -4e-3 250  refr.ind abs.co si.out 5  > /dev/null
  

#	The graphs can be plotted with gnuplot:

#	set par
#	set con
#	set hid
#	splot 's.out' using 1:6:2 w l
#	splot 's.out' using 1:6:3 w l

# here is the  propagation of a non-axial gauss beam:

begin 1e-3 1e-6 100 |  gauss 1.13e-4 2e-4 |\
steps 4e-3 250  refr.ind  abs.co s_s.out 5 > /dev/null

# here we have two non-axial beams:

begin 1e-3 1e-6 100 |  gauss 1.13e-4 2e-4 >foo1
begin 1e-3 1e-6 100 |  gauss 1.13e-4 -2e-4 >foo2
b_mix foo2 < foo1 |\
steps 4e-3 250  refr.ind  abs.co s_ss.out 5 > /dev/null

rm foo*
echo "\n\nYou need a Fortran compiler to run this script!\n\n" 








