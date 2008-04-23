%  LPMultIntensity
%****************************************************************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPMultIntensity multiplies the field with an intensity profile.
%
%	USAGE:
%	Fout = LPMultIntensity(Intens, Fin);
%
%	with:
%	Intens = array of the same dimension as the input field filled with real numbers (doubles) 
%	this can be a bitmap.
%	Use: Intens=double(imread('bitmap','bmp')); in your m-file to read a windows bitmap file 
%	(for example: 'bitmap.bmp') from disk. 
%****************************************************************************************************

function Fout = LPMultIntensity(Intens, Fin)
    msg = ['To change the intensity by a factor K, make sure that you do:', ...
           '  Fin.F .* sqrt(K)'];
    LPMultHelp('LPMultPhase', msg);
    Fout =  0;
endfunction
