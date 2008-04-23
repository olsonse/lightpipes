%LPMultHelp
%****************************************************************************************************
%
%	LightPipes for Octave Optical Toolbox
%
%	Just a little helper to warn about dead LP-components (as compared to
%	the non-gpl matlab version).  This function is used by various dead
%	functions to print a warning.
%
%****************************************************************************************************
function LPMultHelp(func, msg)
    if (!exist('msg') || isempty(msg))
        msg = '';
    endif
    error(['%s:  This function does nothing.  Just do ''Fin.F .* K''\n', ...
           '     to multiply F by an arbitrary field K\n%s'] ...
          , func, msg);
endfunction
