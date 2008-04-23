% LPBeamMix
%*****************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPBeamMix mixes two input fields (superposition).
%
%	USAGE:
%	Fout = LPBeamMix(Fin1,Fin2 [, ratio]);
%
%	with:
%	Fin1, Fin2  = input fields
%       ratio       = ratio of mixture [Default 1]
%                     mixture is ratio*F1.F + 1/(ratio+1e-14)*F2.F
%	Fout        = output field
%*****************************************************
function F = LPBeamMix(F1, F2, ratio)
    if (LPCompatible(F1,F2) != 1)
        error('Incompatible fields cannot be mixed.');
        usage('LPBeamMix');
        return;
    endif

    if (!exist('ratio') || isempty(ratio))
        ratio = 1;
    elseif (ratio > 1 || ratio <0)
        error('ratio must be between 0 and 1');
    endif

    c1 = ratio; c2 = 1/(ratio + 1e-14);

    F.info = F1.info;
    F.F = c1*F1.F + c2*F2.F;
endfunction
