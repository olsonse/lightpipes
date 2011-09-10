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

DEFUN_DLD (LPLensFresnel, args, nargout,
" LPLensFresnel\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPLensFresnel propagates the field in spherical coordinates.\n"
"\n"
"	USAGE:\n"
"	Fout = LPLensFresnel(f,z,Fin);\n"
"\n"
"	with:\n"
"	f = the focal distance of the lens that determines the curvature of the coordinate system\n"
"	z = propagation distance\n"
"	Fin  = input field struct (see LPBegin)\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() != 3 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_map() || !mapIsValidField(octave_stdout, args(2).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;

        print_usage("LPLensFresnel");
        return octave_value();
    }

    double f   = args(0).double_value();
    double z   = args(1).double_value();
    Field * field = mapToField( args(2).map_value() );
    field->lens_fresnel(f, z);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

