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

DEFUN_DLD (LPStrehl, args, nargout,
"LPStrehl\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPStrehl calculates the Strehl ratio of the field.\n"
"\n"
"	USAGE:\n"
"	SR=LPStrehl(Fin);\n"
"\n"
"	with:\n"
"	SR               = calculated Strehl ratio\n"
"	Fin              = the input field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 1 ||
        !args(0).is_map() || !mapIsValidField(octave_stdout, args(0).map_value())
        ) {
        print_usage("LPStrehl");
        return octave_value();
    }

    Field * field = mapToField( args(0).map_value() );
    double strehl = field->get_strehl();
    delete field;       field = NULL;

    return octave_value(strehl);
}

