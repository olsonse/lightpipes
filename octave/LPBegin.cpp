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

DEFUN_DLD (LPBegin, args, nargout,
"LPBegin\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPBegin Creates a new field F.\n"
"\n"
"	USAGE:\n"
"	F = LPBegin( side_length, lambda [, N ] );\n"
"\n"
"	INPUT:\n"
"	side_length   = side length of grid; area of grid is thus side_length^2\n"
"	lambda = wavelength\n"
"	N      = grid dimension [default 256] \n"
"       OUTPUT:\n"
"       F      = octave struct { F , info }\n"
"       F.F    = array of complex field values\n"
"       F.info = N, side_length, lambda from LPBegin input (as well as some internal values)\n"
"\n"
"*****************************************************\n"
) {
    if (args.length() < 2 || args.length() > 3) {
        print_usage("LPBegin");
        return octave_value();
    }

    double side_length   = args(0).double_value();
    double lambda = args(1).double_value();

    int n_grid = 256;
    if (args.length() == 3) n_grid = args(2).int_value();

    Field field(n_grid, side_length, lambda);

    /* Initialize the field to (1,0) */
    {
        int N2 = SQR(field.info.number);
        for ( int j = 0; j < N2; j++ ) {
            field[j] = std::complex<double>(1.,0.);
        }
    }

    return octave_value(fieldToMap(field));
}

