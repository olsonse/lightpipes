% LightPipes for Octave Optical Toolbox
% Calculates the Fraunhofer diffraction of a round hole.
%

m=1;
nm=1e-9*m;
mm=1e-3*m;
cm=1e-2*m;

lambda=1000*nm;
size=10*cm;
N=150;
R=10*mm;
z=1000*m;
f1=200*m;
f2=-200*m;

F=LPBegin(size,lambda,N);
F=LPCircAperture(R,0,0,F);
F=LPLens(f1,0,0,F);
F=LPLensForvard(f2,z,F);
Int=F.F .* conj(F.F);
figure(1);
mymesh(Int);
