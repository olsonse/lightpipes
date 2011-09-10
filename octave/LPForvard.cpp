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

DEFUN_DLD (LPForvard, args, nargout,
"LPForvard\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPForvard propagates the field a distance using the FFT method.\n"
"\n"
"	USAGE:\n"
"\n"
"	Fout = LPForvard(z,Fin);\n"
"	with:\n"
"	z    = propagation distance\n"
"	Fin  = input field struct (see LPBegin)\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() != 2 ||
        !args(0).is_real_scalar() ||
        !args(1).is_map() || !mapIsValidField(octave_stdout, args(1).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;

        print_usage("LPForvard");
        return octave_value();
    }

    double z   = args(0).double_value();
    Field * field = mapToField( args(1).map_value() );
    field->forvard(z);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

