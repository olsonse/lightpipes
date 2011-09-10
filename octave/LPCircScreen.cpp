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

DEFUN_DLD (LPCircScreen, args, nargout,
"LPCircScreen\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPCircScreen screens the field with a circular screen.\n"
"\n"
"	USAGE:\n"
"	Fout = LPCircScreen(R,x0,y0,Fin);\n"
"\n"
"	with:\n"
"	R = the radius of the screen\n"
"	x0, y0 = the position of the center of the screen\n"
"	Fin  = input field struct (see LPBegin)\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() != 4 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_real_scalar() ||
        !args(3).is_map() || !mapIsValidField(octave_stdout, args(3).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() >= 3 && !args(2).is_real_scalar() ) octave_stdout << "arg 2 is invalid" << std::endl;

        print_usage("LPCircScreen");
        return octave_value();
    }

    double r   = args(0).double_value();
    double x0  = args(1).double_value();
    double y0  = args(2).double_value();
    Field * field = mapToField( args(3).map_value() );
    field->circular_screen(r, x0, y0);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

