% LPPhase
%****************************************************************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPPhase calculates the phase of the field.
%
%	USAGE:
%	Phase=LPPhase(Fin);
%
%	with:
%	Fin = the input field
%	Phase = the output phase
%
%****************************************************************************************************
function Fout = LPPhase(Fin)
    msg = ['To obtain the phase, just use arg(Fin.F)'];
    LPMultHelp('LPPhase', msg);
    Fout =  0;
endfunction
