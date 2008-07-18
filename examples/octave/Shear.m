% LightPipes for Octave Optical Toolbox
%Shearing interferometer with aberrated wavefront.

clear;

m=1;
cm=1e-2*m;
mm=1e-3*m;
nm=1e-9*m;

size=4*cm;
lambda=500*nm;
N=256;
R=1*cm;
f=-20*m;
z1=50*cm;
z2=50*cm;
D=3*mm;
D1=1*mm;
Rplate=0.5;
S(1).ab='Spherical aberration with: n=4, m=0, R=10mm, A=10';
S(2).ab='Coma with: n=3, m=-1, R=10mm, A=10';
S(3).ab='Astigmatism with: n=2, m=2, R=7mm, A=10';
S(4).ab='No aberration';
nZ=[4,3,2,0];
mZ=[0,-1,2,0];
RZ=[10*mm,10*mm,7*mm,0];
AZ=[10,10,10,0];

while 1
   k = menu('CHOOSE ABERRATION',...
       		S(1).ab,...
       		S(2).ab,...
       		S(3).ab,...
       		S(4).ab,...
      		'STOP');
   if k == 5
      break;
   end;
   
	F=LPBegin(size,lambda,N);
	F=LPCircAperture(R,0,0,F);
   if k == 4 
      F=LPLens(f,0,0,F);
   else 
      F=LPZernike(nZ(k),mZ(k),RZ(k),AZ(k),F);
   end;
	F=LPForvard(z1,F);
        F1 = F2 = F;
        F1.F *= Rplate;
        F2.F *= (1 - Rplate);
	F2=LPInterpol(size,N,D,D1,0,1,F2);
	F=LPBeamMix(F1,F2);
        Int = F.F .* conj(F.F);

    figure;
 	imshow(Int);
 	title(S(k).ab);
    drawnow;
end
