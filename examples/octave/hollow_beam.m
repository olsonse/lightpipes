% LightPipes for Octave Optical Toolbox

% load the physical constants database
physical
load_LG
propagate_

N               = 512;                  % Number of pixels
side_length     = N* 15.*microns;       % physical size of SLM
lambda          = 780*nm;               % Wavelength
gaussian_size   = 1.00/sqrt(2)*mm;      % 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm;                % LPForvard step size
f               = 20*cm;                % focal length of SLM lens

%%%%%%%% END Param %%%%%%%%%%%%%%


%%%%%%% Create LG order 
[xx yy] = meshgrid(x=3*(-1:2/(N-1):1), y = x);

lg = LG_xy( l=8, p=0, xx, yy, 2/sqrt(l) );
% normalize lg; we only want the phase information
% lg ./= abs(lg);
%mymesh (x,y,arg(lg))
%%%%%%% END Create LG order 




F1 = LPBegin( side_length, lambda, N );
%F2 = LPGaussAperture( gaussian_size,   0, 0, 1, F1);
F2 = F1;
F3 = LPLens( f, 0, 0,   F2);
F3.F .*= lg;
F4 = propagate(f+cm,.5*cm,                F3,0,0,4);
%F4 = LPForvard( 18.0*cm,                          F3);
%F5 = propagate(4.*cm,.5*cm,                F4,0,0,4);
%F5 = propagate(13.*cm, step_size,               F4,0,0,4);

