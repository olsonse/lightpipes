% LightPipes for Octave Optical Toolbox
%
% My replacement for the mesh function of octave.
% The main difference is that mine does not "gset view"
%
% See the help file for mesh()

function mymesh (x, y, z, ifmt)
  switch (nargin)
      case (2)
          fmt = y;
      case (3)
          fmt = z;
      case (4)
          fmt = ifmt;
      otherwise
          fmt = [];
  endswitch

  % change gnuplot if needed. 
  global dont_change_plot_type;
  if (isempty(dont_change_plot_type) || dont_change_plot_type == 0)
      __gnuplot_set__ nosurf
      %__gnuplot_set__ view 0,0,1.9,1
      __gnuplot_set__ view map
      __gnuplot_set__ nokey

      __gnuplot_set__ palette gray
      __gnuplot_set__ pm3d at b
      __gnuplot_set__ noclabel
      __gnuplot_set__ size square
      __gnuplot_set__ noxtics
      __gnuplot_set__ noytics
      __gnuplot_set__ noztics
      __gnuplot_set__ nocolorbox
      dont_change_plot_type = -1; % note already done
  end




% if (ischar(fmt))
  if (isstr(fmt))
    %nargin--;
    nargin = nargin-1;
  else
      fmt = '';
  endif

  if (nargin == 1)
    z = x;
    if (ismatrix (z))
      %gset hidden3d;
      __gnuplot_set__ data style lines;
      %gset surface;
      %gset nocontour;
      __gnuplot_set__ noparametric;
      %gset nologscale;
      %gset view 60, 30, 1, 1
      eval(sprintf('__gnuplot_splot__ (z'') %s',fmt));
    else
      error ("mymesh: argument must be a matrix");
    endif
  elseif (nargin == 3)
    if (isvector (x) && isvector (y) && ismatrix (z))
      xlen = length (x);
      ylen = length (y);
      if (xlen == columns (z) && ylen == rows (z))
        if (rows (y) == 1)
          y = y';
        endif
        len = 3 * xlen;
        zz = zeros (ylen, len);
        k = 1;
        for l = 1:3:len
          zz(:,l)   = x(k) * ones (ylen, 1);
          zz(:,l+1) = y;
          zz(:,l+2) = z(:,k);
          k++;
        endfor
        %gset hidden3d;
        __gnuplot_set__ data style lines;
        %gset surface;
        %gset nocontour;
        %gset nologscale;
        __gnuplot_set__ parametric;
        %gset view 60, 30, 1, 1
        eval(sprintf('__gnuplot_splot__ (zz) %s',fmt));
        __gnuplot_set__ noparametric;
      else
        msg = "mymesh: rows (z) must be the same as length (y) and";
        msg = sprintf ("%s\ncolumns (z) must be the same as length (x)", msg);
        error (msg);
      endif
    elseif (ismatrix (x) && ismatrix (y) && ismatrix (z))
      xlen = columns (z);
      ylen = rows (z);
      if (xlen == columns (x) && xlen == columns (y) &&
        ylen == rows (x) && ylen == rows(y))
        len = 3 * xlen;
        zz = zeros (ylen, len);
        k = 1;
        for l = 1:3:len
          zz(:,l)   = x(:,k);
          zz(:,l+1) = y(:,k);
          zz(:,l+2) = z(:,k);
          k++;
        endfor
        %gset hidden3d;
        __gnuplot_set__ data style lines;
        %gset surface;
        %gset nocontour;
        %gset nologscale;
        __gnuplot_set__ parametric;
        %gset view 60, 30, 1, 1
        eval(sprintf('__gnuplot_splot__ (zz) %s',fmt));
        __gnuplot_set__ noparametric;
      else
        error ("mesh: x, y, and z must have same dimensions");
      endif
    else
      error ("mesh: x and y must be vectors and z must be a matrix");
    endif
  else
    usage ("mesh (z)");
  endif

endfunction
