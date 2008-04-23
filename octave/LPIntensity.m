%LPIntensity
%****************************************************************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	LPIntensity calculates the intensity of the field.
%
%	USAGE:
%	Intens=LPIntensity(flag,Fin);
%
%	with:
%	flag   = 0: no normalization
%	flag   = 1: normalization to one
%	flag   = 2: normalization to 255 for a bitmap picture
%	Fin    = the input field
%	Intens = the calculated intensity of the input field
%
%****************************************************************************************************
function Int = LPIntensity(flag, F)
    switch (flag)
        case (0)
            Int = F.F .* conj(F.F);
        case (1)
            Int = F.F .* conj(F.F);
            Int /= max(max(Int));
        case (2)
            Int = F.F .* conj(F.F);
            Int /= max(max(Int));
            Int *= 255;
    endswitch

    warning(['LPIntensity:  This function doesn''t do too much;\n', ...
             ' just ''F.F .* conj(F.F)'' and then renormalizing that somehow\n', ...
             ' Be brave and just type the matlab code yourself']);
endfunction
