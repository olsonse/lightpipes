% LightPipes for Octave Optical Toolbox

% creates a Laguerre-Gaussian beam (actually a non pure LG mode by projecting
% a Gaussian onto LG(l=1)) and propagates this through a lens, grating, piece
% of glass 1mm thick, to a distance a little further away than the lens focal
% length.

% everything that isn't found here, is found in the LighPipes/Octave example
% directory (download from http://www.umich.edu/~olsonse/software.html).

% load the physical constants database and some other useful stuff
physical
load_LG
propagate_

N               = 512;                  % Number of pixels
grid_dx         = 15.* microns;         % grid cell size
side_length     = N* grid_dx;           % physical size of SLM
lambda          = 780*nm;               % Wavelength
gaussian_size   = 1.00/sqrt(2)*mm;      % 1/e Intensity width of beam; (note, it's not the 1/e^2)
step_size       = 5.*mm;                % LPForvard step size
f               = 20*cm;                % focal length of SLM lens

lf              = 0.64;                 % fraction of f to travel before going through glass.
r               = 0.05;                 % reflection at each surface of glass.
a0              = 45*arc.degrees;       % incident angle of primary beam
n               = 1.55;                 % index of refraction of glass
a               = sin(a0)/n;            % angle of beams in glass (Snell's law)
k               = n*2*pi/lambda;        % k in glass
d               = 1.*mm / cos(a);       % 1-pass distance travelled through glass of secondary beam
dphi            = 2*k*d;                % relative phase delay on the secondary beam
dx              = 2 * d*sin(a);         % displacement of secondary beam parallel to glass.
dx0             = dx * cos(a0);         % displacement of secondary beam transverse to output beam.

P               = 250*mW;               % beam power

%%%%%%%% END Param %%%%%%%%%%%%%%


%%%%%%% Create LG order 
%[xx yy] = meshgrid(x=3*(-1:2/(N-1):1), y = x);
[xx yy] = meshgrid(x=(side_length/2)*(-1:2/(N-1):1), y = x);

lg = LG_xy( l=8, p=0, xx*6/side_length, yy*6/side_length, 2/sqrt(l) );
lg ./= abs(lg);                         % normalize lg; we only want the phase information
%mymesh (x,y,arg(lg))
%%%%%%% END Create LG order 




F    = LPBegin( side_length, lambda, N );
F    = LPGaussAperture( gaussian_size, 0, 0, 1, F);

Pf   = P/sum(sum(F.F.*conj(F.F)));      % power per grid cell [ W ]
Ef   = sqrt(Pf/grid_dx^2);              % sqrt(I) [ sqrt(W/m^2) ]
Ef  /= sqrt(mW/cm^2);                   % put Ef in units [ sqrt(mW/cm^2) ]
F.F *= Ef;                              % put F in units of sqrt(I)

printf('writing SLM phase...\n'); fflush(stdout);
F    = LPLens( f, 0, 0, F);             % apply lens (physical lens)
Fb   = F;                               % store the undeflected beam
F.F.*= lg;                              % apply lg phase
F    = grating(F,32,0);                 % apply a grating phase.
F.F *= sqrt(0.8); F.F += sqrt(0.2)*Fb.F;% mix the beams back together.
printf('finished writing SLM phase.\n'); fflush(stdout);

F    = LPForvard( lf*f, F);

% traverse a piece of glass and interfere primary beam with 1st internally
% reflected+transmitted beam.
Fr    = LPInterpol(side_length,N,dx0,0,0,1, F);
Fr.F *= exp(i * dphi);                  % add phase delay to secondary beam
F.F  *= (1-r); F.F += (r*r*(1-r))*Fr.F; % mix the beams back together.


xlabel('x (mm)'); ylabel('y (mm)'); gset cblabel 'I (mW/cm^2)';
gset colorbox;
F = propagate((1-lf)*f + 2*cm,.5*cm,              F,0,0,1,[],x*1e3,y*1e3);


%  results:  
% Pure LG mode (+grating+lens only), the darkness is down
% to:  ~0.00001 mW/cm^2
%
% Gaussian beam projected onto LG(l=1) mode:
% with neither the non-diffracted fraction of light nor the secondary beam
% from the glass transmission, the darkness at the center of the trap is down
% to   ~0.00185 mW/cm^2

% Adding the secondary beam tranmission from the glass, the darkness rises to:
%      ~0.0300 mW/cm^2

% Adding also the portion of non-diffracted light (assuming the physical
% lens), the darkness is:
%      ~0.022 mW/cm^2 (must have destructively interfered at the center

