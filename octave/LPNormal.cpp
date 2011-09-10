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

DEFUN_DLD (LPNormal, args, nargout,
"LPNormal\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPNormal normalizes the field.\n"
"\n"
"	USAGE:\n"
"	[Fout,NC] = LPNormal(Fin);\n"
"\n"
"	with:\n"
"	NC = calculated normalization coefficient\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 1 ||
        !args(0).is_map() || !mapIsValidField(octave_stdout, args(0).map_value())
        ) {

        print_usage("LPNormal");
        return octave_value();
    }

    Field * field = mapToField( args(0).map_value() );
    double norm_coeff = 0;
    field->normalize(&norm_coeff);
    octave_value_list retval;
    retval.append(fieldToMap(*field));
    retval.append(octave_value(norm_coeff));
    delete field;       field = NULL;

    return retval;
}

