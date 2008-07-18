% LightPipes for Octave Optical Toolbox
%March 1998. F.A. van Goor.
%TwoHoles.m
%Two Holes Interferometer.
clear;

m=1;
nm=1e-9*m;
mm=1e-3*m;
cm=1e-2*m;

lambda=550*nm;
size=5*mm;
N=300;
R=0.12*mm;
d=0.5*mm;
z=5*cm;

F=LPBegin(size,lambda,N);
F1=LPCircAperture(R,d,0,F);
F2=LPCircAperture(R,-d,0,F);
F=LPBeamMix(F1,F2);
clear F1;
clear F2;
figure(1);
xlabel(''); ylabel('');
gset noxtics
gset noytics
for l=1:10
   F=LPForvard(z,F);
   Int=F.F .* conj(F.F);
   subplot(2,5,l);
   mymesh(Int);
   Str=sprintf('z=%4.1f cm',l*z/cm)
   title(Str);
end
clear F;
