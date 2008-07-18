% LightPipes for Octave Optical Toolbox
% Simulates a unstable resonator
%

clear;

m=1; nm=1e-9*m; mm=1e-3*m; cm=1e-2*m;

lambda=308*nm; size=14*mm; N=100; w=5.48*mm;
f1=-10*m; f2=20*m; L=10*m; Isat=1.0; alpha=1e-4; Lgain=1e4;

F=LPBegin(size,lambda,N); F=LPRandomIntensity(F); F=LPRandomPhase(1,F);
for l=1:10
   F=LPRectAperture(w,w,0,0,0,F);   F=LPGain(Isat,alpha,Lgain,F);
   F=LPLensFresnel(f1,L,F);   F=LPGain(Isat,alpha,Lgain,F);
   F=LPLensFresnel(f2,L,F);
   SR(l)=LPStrehl(F);
   F=LPInterpol(size,N,0,0,0,1,F);
   fprintf('Round trip %d Strehl ratio= %f \n',l,SR(l));
   F2=LPRectScreen(w,w,0,0,0,F);
   Int = F2.F .* conj(F2.F);
   %figure(1);   subplot(2,5,l);   imshow(Int);   
   pause(.3);
   mymesh(Int);
end
F2=LPConvert(F2);
Int = F2.F .* conj(F2.F);

figure(2);
__gnuplot_set__ xtics auto;
__gnuplot_set__ ytics auto;
subplot(2,1,1); plot(SR); xlabel('Number of Roundtrips'); ylabel('Strehl ratio');

figure(3); mymesh(Int); %shading interp;
xlabel(''); ylabel('');
colormap(copper); %axis off;
title('Intensity distribution just behind the outcoupler'); %rotate3d on;

%Far-field calculation:
z=1*m; f=40*m;
ff=z*f/(f-z);
F2=LPLens(f,0,0,F2);
F2=LPLensFresnel(ff,z,F2);
F2=LPConvert(F2);
I2= F2.F .* conj(F2.F);

figure(4);
xlabel(''); ylabel('');
mymesh(I2); %shading interp; 
colormap(copper); %axis off;
title('Intensity distribution in the far field'); %rotate3d on;

