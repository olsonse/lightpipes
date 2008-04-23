% LPRandomPhase
%****************************************************************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPRandomPhase substitutes a random phase in the field.
%
%	USAGE:
%
%	Fout=LPRandomPhase(max,Fin);
%
%
%	with:
%	max = the maximum value of the phase.
%	Fin = the input field
%	Fout = the output field
%
%       Note:  use the rand function directly to set and view the seed.
%
%****************************************************************************************************
function Fout = LPRandomPhase(mx,Fin)
    sz = size(Fin.F);
    Fout.F = abs(Fin.F) .* exp(I*mx*rand(sz));
    Fout.info = Fin.info;
endfunction
