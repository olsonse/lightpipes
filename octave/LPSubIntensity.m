%  LPSubIntensity
%****************************************************************************************************
%  THIS FUNCTION IS DEAD IN THIS OCTAVE VERSION.
%
%	LightPipes for Octave Optical Toolbox
%
%	LPSubIntensity substitutes an intensity profile in the field.
%
%	USAGE:
%	Fout = LPSubIntensity(Intensity, Fin);
%
%	with:
%	Intensity = array of the same dimension as the input field filled with real numbers (doubles) 
%	this can be a bitmap.
%	Use: Intensity=double(imread('bitmap','bmp')); in your m-file to read a windows bitmap file 
%	(for example: 'bitmap.bmp') from disk. 
%****************************************************************************************************
function Fout = LPSubIntensity(K, Fin)
    msg = ['To replace the phase, make sure that you do:', ...
           '  K .* exp(i * arg(Fin.F)); note that K can ', ...
           'either be an array or a scalar'];
    LPMultHelp('LPSubIntensity', msg);
    Fout =  0;
endfunction
