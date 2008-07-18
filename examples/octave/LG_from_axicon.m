% LightPipes for Octave Optical Toolbox

% load the physical constants database
physical
load_LG
propagate_

N               = 512;                  % Number of pixels
side_length     = N* 15.*microns;       % physical size of SLM
lambda          = 780*nm;               % Wavelength
gaussian_size   = 1.00/sqrt(2)*mm;      % 1/e Intensity width of beam; (note, it's not the 1/e^2)
axicon_angle    = 175/180*pi;           % included angle of SLM
axicon_n1       = 1.5;                  % index of refraction of the axicon
step_size       = 5.*mm;                % LPForvard step size

%%%%%%%% END Param %%%%%%%%%%%%%%


%%%%%%% Create LG order 
[xx yy] = meshgrid(x=3*(-1:2/(N-1):1), y = x);

lg = LG_xy( l=2, p=0, xx, yy, 2/sqrt(l) );
% normalize lg; we only want the phase information
lg ./= abs(lg);
%mymesh (x,y,arg(lg))
%%%%%%% END Create LG order 




F1 = LPBegin( side_length, lambda,N );
F2 = LPGaussAperture( gaussian_size,   0, 0, 1, F1);
F3 = LPAxicon( axicon_angle, axicon_n1, 0, 0,   F2);
%F4 = propagate(13.*cm, step_size,               F3,0,0,4);
F4 = LPForvard( 7.5*cm,                          F3);
F5 = LPAxicon( axicon_angle, axicon_n1, 0, 0,   F4);
        printf('moving towards slm\n'); fflush(stdout);
%F6 = LPForvard( 40*cm,                          F5);
F6 = propagate(40.*cm, 10*cm,               F5,0,0,4);
        printf('applying slm\n');       fflush(stdout);
F6.F .*= lg;
F7 = propagate(50.*cm,2*cm,                F6,0,0,4);

