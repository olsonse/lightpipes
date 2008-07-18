// -*- c++ -*-
// $Id$
/*
 * Copyright 2005 Spencer Olson
 *
 * $Log$
 *
 */


#include "oct-Field.h"
#include <octave/oct.h>
#include <octave/pager.h>

DEFUN_DLD (LPInterpol, args, nargout,
"LPInterpol\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPInterpol interpolates the field to a new grid.\n"
"\n"
"	USAGE:\n"
"	Fout=LPInterpol(new_side_length,new_N,x0,y0,angle,magnification,Fin);\n"
"\n"
"	with:\n"
"	new_side_length  = new grid size after interpolation\n"
"	new_N            = new grid dimension after interpolation (grid of NxN elements)\n"
"	x0, y0           = displacement from the center\n"
"	angle            = rotation of the field\n"
"	magnification    = magnification of the field\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 7 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_real_scalar() ||
        !args(3).is_real_scalar() ||
        !args(4).is_real_scalar() ||
        !args(5).is_real_scalar() ||
        !args(6).is_map() || !mapIsValidField(octave_stdout, args(6).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() >= 3 && !args(2).is_real_scalar() ) octave_stdout << "arg 2 is invalid" << std::endl;
        if ( args.length() >= 4 && !args(3).is_real_scalar() ) octave_stdout << "arg 3 is invalid" << std::endl;
        if ( args.length() >= 5 && !args(4).is_real_scalar() ) octave_stdout << "arg 4 is invalid" << std::endl;
        if ( args.length() >= 6 && !args(5).is_real_scalar() ) octave_stdout << "arg 5 is invalid" << std::endl;

        print_usage("LPInterpol");
        return octave_value();
    }

    double L   = args(0).double_value();
    int    N   = args(1).int_value();
    double x0  = args(2).double_value();
    double y0  = args(3).double_value();
    double a   = args(4).double_value();
    double M   = args(5).double_value();

    if (M <=0 || N <=0 || L <= 0) {
        octave_stdout << "new_side_length, new_number, "
                         "and magnification must all be > 0 !"
                      << std::endl;
        print_usage("LPInterpol");
        return octave_value();
    }

    Field * field = mapToField( args(6).map_value() );
    field->interpolate(L, N, x0, y0, a, M);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

