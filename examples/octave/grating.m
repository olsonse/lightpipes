% LightPipes for Octave Optical Toolbox
%%%%%%%%%%%%%%%%%%%% BEGIN GRATING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = grating(F,px,py)
% F = grating(F,px,py)
%   applies a grating with px,py periods in x,y directions respectively.
    Gx = make_grating(F.info.number,px);
    Gy = make_grating(F.info.number,py);

    for k = 1:F.info.number
        F.F(:,k) *= Gx(k);
        F.F(k,:) *= Gy(k);
    end
endfunction


function G = make_grating(N,p)
% G = grating(p)
%   returns the field of a grating function with 'p' periods over the 512
%   pixels.
    G = (0:(N-1)) * ( 2*pi * p / N );
    G = exp(i*G);
endfunction


%%%%%%%%%%%%%%%%%%%%%% END GRATING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
