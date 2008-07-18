% LightPipes for Octave Optical Toolbox

printf(['Various useful functions are now defined:\n', ...
        '\t LG, LG_xy, fact, PHI, LaguerreL\n']  ...
);

function lg = LG(l,p,rho,phi,omega0)
% function lg = LG(l,p,rho,phi,omega0)
%       returns the field of a laguerre gaussian beam of order (l,p).  Note
%       that the result is unitless (i.e. scalable in position space).
%       where:
%               l   = azimuthal order
%               p   = radial order
%               rho = radial position of evaluation of LG mode
%               phi = azimuthal position of evaluation of LG mode
%               omega0 = 
%       Note that the electric field amplitude is given by
%           (sqrt(P_{0})/omega0) * LG(l,p,rho,phi,omega0)
%       where P_{0} is the beam power.
    lg = ...
    sqrt(2 * fact(p) / (pi * fact(p+abs(l))))   ...
    .* (sqrt(2) * rho / omega0).^(abs(l))        ...
    .* LaguerreL( abs(l), p, 2*rho.^2 / omega0^2 )         ...
    .* exp( -rho.^2/omega0^2 + i*l*phi);
endfunction

function angle = PHI(x,y)
% function angle = PHI(x,y)
%       return the positive (CCW) angle from the positive x axis.  This
%       function correctly determines the angle independent of quadrant.
%       where:
%               x = x position
%               y = y position
    angle = acos(x./sqrt(x.^2 + y.^2));
    k = find(y < 0);
    angle(k) = 2*pi - angle(k);
endfunction;


function lg = LG_xy(l,p,xx,yy,omega0)
% function lg = LG_xy(l,p,xx,yy,omega0)
%       where:
%               l   = azimuthal order
%               p   = radial order
%               xx  = x position of evaluation of LG mode
%               yy  = y position of evaluation of LG mode
%               omega0 = 
%       Note that this function is defined in terms of LG(...) and PHI(x,y)
    lg = LG(l,p,sqrt(xx.^2+yy.^2), PHI(xx,yy), omega0);
endfunction

function retval = fact (n)
% function retval = fact (n)
%       returns n! where n is an integer
    n = floor(n);
    if (n > 0)
%       retval = n * fact (n-1);

        retval = 1;
        while (n >= 1)
            retval *= n;
            --n;
        end
    else
        retval = 1;
    endif
endfunction

function genL = LaguerreL(l,p,rho)
% function genL = LaguerreL(l,p,rho)
%       returns the general Laguerre function L(l,p,rho)
%       where:
%               l = azimuthal order
%               p = radial order
%               rho = radial position
    genL = 0;
    for m = 0:p
        genL += (-1)^m / (fact(p-m) * fact(abs(l)+m) * fact(m)) * rho.^m;
    end

    genL *= fact(abs(l)+p);
endfunction


