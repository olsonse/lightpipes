% LightPipes for Octave Optical Toolbox
% calculation of the intensity near the
% focus of a lens with LPSteps.
% F.A. van Goor

clear;

m=1;
nm=1e-9*m;
mm=1e-3*m;
cm=1e-2*m;

lambda=632.8*nm;
size=4*mm;
N=100;
R=1.5*mm;
dz=10*mm;
f=50*cm;
n=(1+0.1*i)*ones(N,N);

F=LPBegin(size,lambda,N);
F=LPCircAperture(R,0,0,F);
F=LPLens(f,0,0,F);
for l=1:100
   F=LPSteps(dz,1,n,F);
   Int = F.F .* conj(F.F);
   for k=1:N
      Icross(l,k)=Int(N/2,k);
   end
end
figure;
mymesh(Icross);
