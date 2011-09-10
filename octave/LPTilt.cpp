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

DEFUN_DLD (LPTilt, args, nargout,
"LPTilt\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPTilt introduces a tilt in the field.\n"
"\n"
"	USAGE:\n"
"	Fout=LPTilt(tx,ty,Fin);\n"
"\n"
"	with:\n"
"	tx, ty = tilt angles in horizontal and vertical directions in radians\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 3 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_map() || !mapIsValidField(octave_stdout, args(2).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;

        print_usage("LPTilt");
        return octave_value();
    }

    double tx  = args(0).double_value();
    double ty  = args(1).double_value();
    Field * field = mapToField( args(2).map_value() );
    field->tilt(tx, ty);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

