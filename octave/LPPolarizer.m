% LPPolarizer
%****************************************************************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPPolarizer outputs a linear polarized field.
%
%	USAGE:
%	Fout=LPPolarizer(Phi,Fx,Fy);
%
%	with:
%	Phi=polarizer angle in radians (with respect to the positive x-axis)
%	Fx, Fy = input fields of horizontal and vertical components
%	Fout   = the output field
%
%****************************************************************************************************
function F = LPPolarizer(Phi, Fx, Fy)
    F = LPBeamMix(Fx, Fy, cos(Phi));
endfunction
