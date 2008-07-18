% LightPipes for Octave Optical Toolbox

% This file is for taking a look at the optical field for doing evaporative
% cooling by manufacturing a major leak at the top of the trap.  

% This version looks at the result of placing the knife edge between the two
% relaying lenses (that relay the SLM plane to fslm away from the MOT).

% load the physical constants database
physical
load_LG
propagate_

N               = 512;                  % Number of pixels
side_length     = N* 15.*microns;       % physical size of SLM
lambda          = 780*nm;               % Wavelength
gaussian_size   = 0.40/sqrt(2)*mm;      % 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm;                % LPForvard step size
fslm            = 15*cm;                % focal length of SLM lens
fslmR           = 30*cm;                % focal length of first relay (to move SLM plane)
fr1             = 20*cm;                % focal length of relay:2 lens:[1,4] (around MOT)
fr2             =  8*cm;                % focal length of relay:2 lens:[2,3] (around MOT)

%%%%%%%% END Param %%%%%%%%%%%%%%


%%%%%%% Create LG order 
[xx yy] = meshgrid(x=3*(-1:2/(N-1):1), y = x);

lg = LG_xy( l=8, p=0, xx, yy, 2/sqrt(l) );
% normalize lg; we only want the phase information
lg ./= abs(lg);
%mymesh (x,y,arg(lg))
%%%%%%% END Create LG order 



if (0)
%%%%%  Just the hollow beam relayed. 
F1 = LPBegin(side_length, lambda, N );                  % beginning field
F2 = LPGaussAperture(gaussian_size, 0, 0, 1, F1);
F2.F .*= lg;                                            % SLM LG phase
F2 = LPLens(    fslm,      0, 0, F2);                    % 'SLM' lens
F3 = LPForvard( fslmR,           F2);
F4 = LPLens(    fslmR,     0, 0, F3);                    % relay:1 lens:1
F5 = LPForvard(2*fslmR,          F4);
F6 = LPLens(    fslmR,     0, 0, F5);                    % relay:1 lens:2
F7 = LPForvard( fslmR+fslm,      F6);                    % MOT (1st pass)
    mymesh(F7.F .* conj(F7.F));
    fprintf (stderr, 'press enter to continue'); pause;
F8 = LPForvard( fr1,             F7);
F9 = LPLens(    fr1,       0, 0, F8);                    % relay:2 lens:1
F10= LPForvard( fr1+fr2,         F9);
F11= LPLens(    fr2,       0, 0, F10);                   % relay:2 lens:2
F12= LPForvard( 2*fr2,           F11);
F13= LPLens(    fr2,       0, 0, F12);                   % relay:2 lens:3
F14= LPForvard( fr1+fr2,         F13);
F15= LPLens(    fr1,       0, 0, F14);                   % relay:2 lens:1
F16= LPForvard( fr1,             F15);                   % MOT (2nd pass)
%F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4);
    mymesh(F16.F .* conj(F16.F));
end



if (0)
%%%%%  With the knife edge between the SLM (at SLM focus) relay:1,lens:1.
F1 = LPBegin(side_length, lambda, N );                  % beginning field
F2 = LPGaussAperture(gaussian_size, 0, 0, 1, F1);
F2.F .*= lg;                                            % SLM LG phase
F2 = LPLens(    fslm,      0, 0, F2);                    % 'SLM' lens
F3 = LPForvard( fslm,            F2);
F3 = LPRectScreen(50,100,25+100*microns,0,0,F3);        % knife edge
F3 = LPForvard( fslmR-fslm,      F3);
F4 = LPLens(    fslmR,     0, 0, F3);                    % relay:1 lens:1
F5 = LPForvard( 2*fslmR,         F4);
F6 = LPLens(    fslmR,     0, 0, F5);                    % relay:1 lens:2
F7 = LPForvard( fslmR+fslm,      F6);                    % MOT (1st pass)
    mymesh(F7.F .* conj(F7.F));
    fprintf (stderr, 'press enter to continue'); pause;
F8 = LPForvard( fr1,             F7);
F9 = LPLens(    fr1,       0, 0, F8);                    % relay:2 lens:1
F10= LPForvard( fr1+fr2,         F9);
F11= LPLens(    fr2,       0, 0, F10);                   % relay:2 lens:2
F12= LPForvard( 2*fr2,           F11);
F13= LPLens(    fr2,       0, 0, F12);                   % relay:2 lens:3
F14= LPForvard( fr1+fr2,         F13);
F15= LPLens(    fr1,       0, 0, F14);                   % relay:2 lens:1
F16= LPForvard( fr1,             F15);                   % MOT (2nd pass)
%F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4);
    mymesh(F16.F .* conj(F16.F));
end

if (1)
%%%%%  With the knife edge between the SLM (at SLM focus) relay:1,lens:1.
%%%%%  In this case, we relay the knife edge to the MOT and not the SLM plane.
F1 = LPBegin(side_length, lambda, N );                  % beginning field
F2 = LPGaussAperture(gaussian_size, 0, 0, 1, F1);
F2.F .*= lg;                                            % SLM LG phase
F2 = LPLens(    fslm,      0, 0, F2);                    % 'SLM' lens
F3 = LPForvard( fslm,            F2);
F3 = LPRectScreen(50,100,25+100*microns,0,0,F3);        % knife edge
F3 = LPForvard( fslmR,           F3);
F4 = LPLens(    fslmR,     0, 0, F3);                    % relay:1 lens:1
F5 = LPForvard( 2*fslmR,         F4);
F6 = LPLens(    fslmR,     0, 0, F5);                    % relay:1 lens:2
F7 = LPForvard( fslmR,           F6);                    % MOT (1st pass)
    mymesh(F7.F .* conj(F7.F));
    fprintf (stderr, 'press enter to continue'); pause;
F8 = LPForvard( fr1,             F7);
F9 = LPLens(    fr1,       0, 0, F8);                    % relay:2 lens:1
F10= LPForvard( fr1+fr2,         F9);
F11= LPLens(    fr2,       0, 0, F10);                   % relay:2 lens:2
F12= LPForvard( 2*fr2,           F11);
F13= LPLens(    fr2,       0, 0, F12);                   % relay:2 lens:3
F14= LPForvard( fr1+fr2,         F13);
F15= LPLens(    fr1,       0, 0, F14);                   % relay:2 lens:1
F16= LPForvard( fr1,             F15);                   % MOT (2nd pass)
%F16= propagate( fslm, 0.1*cm, F16, 0, 0, 4);
    mymesh(F16.F .* conj(F16.F));
end

