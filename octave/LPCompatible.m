% LPCompatible
%*****************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPCompatible determines if two fields can be mixed.
%
%	USAGE:
%	compatible = LPCompatible(Fin1,Fin2);
%
%	with:
%	Fin1, Fin2  = input fields
%	compatible  = output field
%*****************************************************
function compatible = LPCompatible(F1,F2)
    if (F1.info.number == F2.info.number && ...
        F1.info.lambda == F2.info.lambda && ...
        F1.info.side_length == F2.info.side_length && ...
        F1.info.sph_coords_factor == F2.info.sph_coords_factor)
        compatible = 1;
    else
        compatible = 0;
    endif
endfunction
