% LightPipes for Octave Optical Toolbox
%March 1998. F.A. van Goor.
%OneHole.m
%One Hole Diffraction.

clear;

m=1;
nm=1e-9*m;
mm=1e-3*m;
cm=1e-2*m;

lambda=550*nm;
size=5*mm;
N=100;
R=1*mm;
z=25*cm;

F = LPBegin(size,lambda,N);
F = LPCircAperture(R,0,0,F);
F = LPRectScreen(R*2,R/8,0,0,-45,F);
F = LPCircScreen(R/4,0,0,F);
F = LPFresnel(z,F);
Int = F.F .* conj(F.F);
%imshow(Int);
title('Intensity Distribution in the somewhat far-field');
mymesh(Int);
