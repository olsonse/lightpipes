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

DEFUN_DLD (LPAxicon, args, nargout,
"LPAxicon\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPAxicon simulates an axicon.\n"
"\n"
"	USAGE:\n"
"	Fout = LPAxicon(phi,n1,x0,y0,Fin);\n"
"\n"
"	with:\n"
"	phi = top angle of the axicon in radians\n"
"	n1  = complex refractive index of the axicon material\n"
"	x0, y0 = position of the center\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 5 ||
        !args(0).is_real_scalar() ||
        !args(1).is_scalar_type() ||
        !args(2).is_real_scalar() ||
        !args(3).is_real_scalar() ||
        !args(4).is_map() || !mapIsValidField(octave_stdout, args(4).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_scalar_type() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() >= 3 && !args(2).is_real_scalar() ) octave_stdout << "arg 2 is invalid" << std::endl;
        if ( args.length() >= 4 && !args(3).is_real_scalar() ) octave_stdout << "arg 3 is invalid" << std::endl;

        print_usage();
        return octave_value();
    }

    double phi   = args(0).double_value();
    std::complex<double>
           n1    = args(1).complex_value();
    double x0    = args(2).double_value();
    double y0    = args(3).double_value();
    Field * field = mapToField( args(4).map_value() );
    field->axicon(phi, n1, x0, y0);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

