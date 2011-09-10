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

DEFUN_DLD (LPConvert, args, nargout,
"LPConvert\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPConvert converts from spherical coordinates to normal coordinates.\n"
"\n"
"	USAGE:\n"
"	Fout = LPConvert(Fin);\n"
"\n"
"	with:\n"
"	Fin  = input field struct (see LPBegin)\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() != 1 ||
        !args(0).is_map() || !mapIsValidField(octave_stdout, args(0).map_value())
        ) {

        print_usage("LPConvert");
        return octave_value();
    }

    Field * field = mapToField( args(0).map_value() );
    field->spherical_to_normal_coords();
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

