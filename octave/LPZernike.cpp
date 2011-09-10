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

DEFUN_DLD (LPZernike, args, nargout,
"LPZernike\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPZernike introduces arbitrary Zernike aberration into the field.\n"
"\n"
"	USAGE:\n"
"	Fout=LPZernike(n,m,R,A,Fin);\n"
"\n"
"	with:\n"
"	n, m             = the integer orders,\n"
"	        (See Born and Wolf, 6th edition p.465, Pergamon 1993)\n"
"	R                = the radius at which the phase amplitude reaches A  radians\n"
"	A                = phase amplitude\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 5 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_real_scalar() ||
        !args(3).is_real_scalar() ||
        !args(4).is_map() || !mapIsValidField(octave_stdout, args(4).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() >= 3 && !args(2).is_real_scalar() ) octave_stdout << "arg 2 is invalid" << std::endl;
        if ( args.length() >= 4 && !args(3).is_real_scalar() ) octave_stdout << "arg 3 is invalid" << std::endl;

        print_usage("LPZernike");
        return octave_value();
    }

    int n       = args(0).int_value();
    int m       = args(1).int_value();
    double R    = args(2).double_value();
    double A    = args(3).double_value();
    Field * field = mapToField( args(4).map_value() );
    field->zernike(n, m, R, A);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

