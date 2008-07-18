% LightPipes for Octave Optical Toolbox
% Twyman Green interferometer.
% The Twyman green interferometer consists of a beamsplitter
% and two mirrors. One of the mirrors
% is perfect, the other can be investigated for aberrations. 
%
% F.A. van Goor, March 1998.

clear;

m=1;
cm=1e-2*m;
mm=1e-3*m;
nm=1e-9*m;

size=30*mm;
lambda=500*nm;
R=5*mm;
N=250;
z1=50*cm;
z2=40*cm;
z3=40*cm;
z4=100*cm;
RBS=0.3;
nz=3;
mz=1;
Rz=0.005;
Az=25;

F=LPBegin(size,lambda,N); F=LPCircAperture(R,0,0,F);
F=LPForvard(z1,F);
F1 = F2 = F; F1.F *= RBS; F2.F *= 1-RBS;
F1=LPForvard(2*z2,F1); F2=LPForvard(z3,F2);
F2=LPZernike(nz,mz,Rz,Az,F2);
F2=LPForvard(z3,F2);
F1.F *= RBS; F2.F *= 1-RBS;
F=LPBeamMix(F1,F2);
F=LPForvard(z4,F); F=LPInterpol(0.012,250,0,0,0,1,F);
Int=F.F .* conj(F.F);

figure(1);
mymesh(Int)
title('Twyman Green interferometer with Zernike aberration');
