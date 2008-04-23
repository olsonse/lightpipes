%LPIntAttenuator
%*****************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPIntAttenuator attenuates the field's intensity by a factor R.
%
%	USAGE:
%	Fout=LPIntAttenuator(R,Fin);
%
%	with:
%	R  = attenation
%	Fin  = input field
%	Fout = output field
%*****************************************************
function Fout = LPIntAttenuator(R, Fin)
    msg = ['To change the intensity by a factor K, make sure that you do:', ...
           '  Fin.F .* sqrt(K)'];
    LPMultHelp('LPIntAttenuator', msg);
    Fout =  0;
endfunction
