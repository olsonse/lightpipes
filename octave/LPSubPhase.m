% LPSubPhase
%****************************************************************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPSubPhase substitutes a phase profile in the field.
%
%	USAGE:
%	Fout = LPSubPhase(Phase, Fin);
%
%	with:
%	Phase = array of the same dimension as the input field filled with real numbers (doubles) 
%	this can be a bitmap.
%	Use: Phase=double(imread('bitmap','bmp')); in your m-file to read a windows bitmap file 
%	(for example: 'bitmap.bmp') from disk. 
%****************************************************************************************************
function Fout = LPSubPhase(Phase, Fin)
    msg = ['To replace the phase, make sure that you do:', ...
           '  abs(Fin.F) .* exp(i * K); note that K can ', ...
           'either be an array or a scalar'];
    LPMultHelp('LPSubPhase', msg);
    Fout =  0;
endfunction

