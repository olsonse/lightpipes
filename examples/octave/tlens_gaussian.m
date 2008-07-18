% LightPipes for Octave Optical Toolbox

% load physical units and constants
physical  

% load some utility functions for watching propagation
propagate_

% pixel size on slm
DX=15*microns;
% wavelength
lambda=780*nm;

% lens applied directly on camera.
f_torus=40*cm;

% The radius of the center of the toroidal lens.
R0=128*DX;

R1 = .19*512*DX; % central radius of phase ring.
w1 = .23*512*DX; % width of phase ring.

%w = 0.12*40*mm/sqrt(2); % 40mm colimating lens
w = 1.71*mm/sqrt(2); % smaller beam

function phi = tlens(x,y,R0,f,lambda);
% function phi = tlens(x,y,R0,f,lambda);
%
%  This function calculates the phase to use for a
%   toroidal lens + quadratic lens.
    phi = 2*pi + rem(-pi*(sqrt(x.^2+y.^2) - R0).^2 / (f * lambda), 2*pi);
endfunction

function phi = phase_ring(x, y, R, w)
% function phi = phase_ring(x, y, R, w)
%
%  This function calculates the phase ring that is to be added to the phase of
%  the toroidal lens at the SLM.
    phi = 0 * x;
    r2 = (x.^2 + y.^2);
    Ri_2 = (R-0.5*w)^2.0;
    Rf_2 = (R+0.5*w)^2.0;
    phi(find( r2 >= Ri_2 & r2 < Rf_2 )) = pi;
endfunction

function plot_phase(phase, div)
% just a helper function to plot the phase and subsample first to get gnuplot
% to quit taking forever.
    gset view map
    gset nosurf
    gset nohidden3d
    gset pm3d at b

    gsplot phase(1:div:rows(phase),1:div:columns(phase));
endfunction


[xx yy] = meshgrid( x=(-255.5:1:255.5)*DX, y=x);
torus_phase = tlens (xx,yy,R0,f_torus,lambda) + phase_ring(xx,yy,R1,w1);

%plot_phase(torus_phase,4);
input('press enter to continue');



% propagate the field from the SLM and take a look
F1 = LPBegin(512*DX, lambda, 512);              % initial field
F2 = LPGaussAperture(w, 0, 0, 1,        F1);    % gaussian input
F3 = LPCircAperture(256*DX,0,0,         F2);    % Iris in front of SLM
F4 = F3; F4.F .*= exp(i * torus_phase );        % apply SLM phase
F5 = LPForvard( f_torus,                F4 );
%F5 = propagate(f_torus, f_torus/25,     F4, 0, 0, 4, 256);

plot_phase(F5.F .* conj(F5.F),4);
