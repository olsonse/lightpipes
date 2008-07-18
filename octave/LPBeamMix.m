% LPBeamMix
%*****************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPBeamMix mixes two input fields (superposition).
%
%	USAGE:
%	Fout = LPBeamMix(Fin1,Fin2);
%
%	with:
%	Fin1, Fin2  = input fields
%	Fout        = output field
%*****************************************************
function F = LPBeamMix(F1, F2)
    if (LPCompatible(F1,F2) != 1)
        error('Incompatible fields cannot be mixed.');
        usage('LPBeamMix');
        return;
    endif

    F = F1;
    F.F += F2.F;
endfunction
