% LPMultPhase
%****************************************************************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPMultPhase multiplies the field with an phase profile.
%
%	USAGE:
%	Fout = LPMultPhase(Phase, Fin);
%
%	with:
%	Phase = array of the same dimension as the input field filled with real numbers (doubles) 
%	this can be a bitmap.
%	Use: Phase=double(imread('bitmap','bmp')); in your m-file to read a windows bitmap file 
%	(for example: 'bitmap.bmp') from disk. 
%****************************************************************************************************

function Fout = LPMultPhase(Phase, Fin)
    msg = ['To change the phase by a factor K, make sure that you do:', ...
           '  Fin.F .* exp(i * K)'];
    LPMultHelp('LPMultPhase', msg);
    Fout =  0;
endfunction
