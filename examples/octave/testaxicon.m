% LightPipes for Octave Optical Toolbox

gset nosurf
gset nohidden3d
gset pm3d at b
gset view map
axis([0 256 0 256]);

phi = 178 * (pi/180);


F0  = LPBegin(512*15e-6, 780e-9, 256);
F1  = LPGaussAperture(.75e-3, 0, 0, 1,          F0); 
F2  = LPAxicon (phi, 1.5,   0e-6, 0,            F1);
F3  = LPForvard( 25e-2,                         F2);
%   gsplot F3.F.*conj(F3.F);
%   input ('press enter to continue');
F4  = LPAxicon (phi, 1.5,   0e-6, 0,            F3);
F5  = LPForvard( 100e-2,                        F4);
    gsplot F5.F.*conj(F5.F);
%   input ('press enter to continue');

%F6  = LPLens( -20e-2, 0, 0,                     F5);
%F7  = LPForvard(  3e-2,                         F6);
%    gsplot F7.F.*conj(F7.F);
%
