% LightPipes for Octave Optical Toolbox
0;

printf('functions defined : propagate, propagate2fields\n'); fflush(stdout);

function F = propagate(z, step_size, F, plot_arg, spherical, ds, plotline, x, y)
% function F = propagate(z, step_size, F, plot_arg, spherical, ds, plotline,x,y)
%  where
%       z               : total distance to propagate
%       step_size       : distance to propagate between plots
%       F               : the input field (see LPBegin)
%       plot_arg        : 0 : plot intensity, 1 : plot phase
%       spherical       : attempt the LPLensForvard
%       ds              : subsampling for plotting
%       plotline        : -1 : plot surface, 1-columns(F) : plot column of F
%       x               : x coords used in mymesh
%       y               : y coords used in mymesh
    if (!exist('plot_arg') || isempty(plot_arg))
        plot_arg = 0;
    endif

    if (!exist('spherical') || isempty(spherical))
        spherical = 0;
    endif

    if (!exist('ds') || isempty(ds))
        ds = 2;
    endif

    if (!exist('plotline') || isempty(plotline))
        plotline = -1;
    endif

    if (!exist('x') || isempty(x))
        x = 1:rows(F.F);;
    endif

    if (!exist('y') || isempty(y))
        y = 1:columns(F.F);
    endif

    for j=0:step_size:z;
        printf('%g, ', j); fflush(stdout);

        if (spherical == 0)
            F=LPForvard(step_size,F);
        else
            F=LPLensForvard(spherical,step_size,F);
            F=LPConvert(F);
        endif

        if (plotline < 0)
            val = F.F(1:ds:F.info.number,1:ds:F.info.number);
            xx = x(1:ds:F.info.number);
            yy = y(1:ds:F.info.number);

            if (plot_arg == 0)
                val = abs(val.*conj(val));
            else
                val = arg(val);
            endif
            mymesh(xx,yy,val);
        else
            fl=F.F(:,plotline);
            if (plot_arg == 0)
                fl.*=conj(fl);
            else 
                fl=arg(fl);
            endif
            plot(fl)
        end
    end;
    printf('\n');
endfunction


function [F F2] = propagate2fields(z, step_size, F, F2, ratio, plot_arg, spherical, ds, prop1, prop2)
% function F = propagate2fields(z, step_size, F, plot_arg, spherical, ds, plotline)
%       This function is for copropagating two fields, e.g. two polarizations
%  where
%       z               : total distance to propagate
%       step_size       : distance to propagate between plots
%       F               : the input field (see LPBegin)
%       F2              : a second input field (see LPBegin)
%       ratio           : ratio*|F|^2 + (1-ratio)*|F2|^2 is plotted
%       plot_arg        : 0 : plot intensity, 1 : plot phase
%       spherical       : attempt the LPLensForvard
%       ds              : subsampling for plotting
%       plotline        : 0 : plot surface, 1 : plot line
    if (!exist('plot_arg') || isempty(plot_arg))
        plot_arg = 0;
    endif

    if (!exist('spherical') || isempty(spherical))
        spherical = 0;
    endif

    if (!exist('ds') || isempty(ds))
        ds = 2;
    endif

    if (!exist('ratio') || isempty(ratio))
        ratio = 1;
    endif

    if (!exist('prop1') || isempty(prop1))
        prop1 = 1;
    endif

    if (!exist('prop2') || isempty(prop2))
        prop2 = 1;
    endif

    for j=0:step_size:z;
        printf('%g, ', j); fflush(stdout);

        if (spherical == 0)
            if (prop1 == 1)
                F=LPForvard(step_size,F);
            endif
            if (prop2 == 1)
                F2=LPForvard(step_size,F2);
            endif
        else
            if (prop1 == 1)
                F=LPLensForvard(spherical,step_size,F);
                F=LPConvert(F);
            endif
            if (prop2 == 1)
                F2=LPLensForvard(spherical,step_size,F2);
                F2=LPConvert(F2);
            endif
        endif

        val1 = F.F(1:ds:F.info.number,1:ds:F.info.number);
        val2 = F2.F(1:ds:F2.info.number,1:ds:F2.info.number);

        if (plot_arg == 0)
            val = ratio*val1.*conj(val1) + (1-ratio)*val2.*conj(val2);
            val = abs(val);
        else
            % it doesn't really make sense to plot the arg of both, so we'll
            % just do val1
            val = arg(val1);
        endif
        mymesh(val);
    end;
    printf('\n');
endfunction
