% LightPipes for Octave Optical Toolbox

load_LG;
propagate_;

lg = LG_xy( l=6, p=0, xx, yy, 2/sqrt(l) );
lg ./= abs(lg);



f1  = F2; f1.F .*= lg;
f1  = LPLens(-150*cm,0,0,                       f1);
fg1 = F2;
%   f2 = propagate(100.*cm,2*cm,                    f1);
%   fg2 = propagate(100.*cm,2*cm,                   F2);
[f2 fg2] = propagate2fields(100*cm,2*cm,              f1,fg1,0.9,0,0,2,1,1);
f3  = LPLens(-20*cm,0,0,                        f2);
fg3 = LPLens(-20*cm,0,0,                        fg2);
%   f4 = propagate(20.*cm,step_size,                f3);
%   fg4 = propagate(20.*cm,step_size,               fg3);
[f4 fg4] = propagate2fields(20*cm,step_size,          f3,fg3,.9,0,0,2,1,1);


[f5 fg5] = propagate2fields(1.375.*cm,.25*step_size,  f4,fg4,.9,0,0,2,1,0);

