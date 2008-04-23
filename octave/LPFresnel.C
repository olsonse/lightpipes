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

DEFUN_DLD (LPFresnel, args, nargout,
"LPFresnel\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPFresnel propagates the field a distance using the convolution method.\n"
"\n"
"	USAGE:\n"
"	Fout = LPFresnel(z,Fin);\n"
"\n"
"	with:\n"
"	z    = propagation distance\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 2 ||
        !args(0).is_real_scalar() ||
        !args(1).is_map() || !mapIsValidField(octave_stdout, args(1).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;

        print_usage("LPFresnel");
        return octave_value();
    }

    double z   = args(0).double_value();
    Field * field = mapToField( args(1).map_value() );
    field->fresnel(z);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

