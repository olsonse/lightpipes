% LPRandomIntensity
%****************************************************************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPRandomIntensity substitutes a random intensity in the field.
%       This random intensity is normalized to 1.  Use Fout *= 1e-3 to
%       normalize to 1 milliwatt for example.
%
%	USAGE:
%	Fout=LPRandomIntensity(Fin);
%
%	with:
%	seed = an arbitrary number to initiate the random number generator
%	Fin = the input field
%	Fout = the output field
%
%       Note:  use the rand function directly to set and view the seed.
%
%****************************************************************************************************

function Fout = LPRandomIntensity(Fin)
    sz = size(Fin.F);
    Fout.F = Fin.F ./ abs(Fin.F) .* rand(sz);
    Fout.info = Fin.info;
endfunction
