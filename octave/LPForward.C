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

DEFUN_DLD (LPForward, args, nargout,
"LPForward\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPForward propagates the field a distance using the direct integration method.\n"
"\n"
"	USAGE:\n"
"	Fout = LPForward(z, NewSize, NewNumber, Fin);\n"
"\n"
"	with:\n"
"	z    = propagation distance\n"
"	NewSize = a new grid size\n"
"	NewNumber = a new grid dimension\n"
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

        print_usage("LPForward");
        return octave_value();
    }

    double z   = args(0).double_value();
    double sz  = args(1).double_value();
    int    N   = args(2).int_value();
    Field * field = mapToField( args(3).map_value() );
    field->forward(z, sz, N);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

