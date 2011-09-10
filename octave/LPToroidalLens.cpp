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

DEFUN_DLD (LPToroidalLens, args, nargout,
"LPToroidalLens\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPToroidalLens propagates the field through a toroidal lens.\n"
"\n"
"	USAGE:\n"
"	Fout=LPToroidalLens(R,f,x0,y0,Fin);\n"
"\n"
"	with:\n"
"	R                = toroidal radius\n"
"	f                = focal length of the lens\n"
"	x0, y0           = displacement from the center\n"
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

        print_usage("LPToroidalLens");
        return octave_value();
    }

    double R   = args(0).double_value();
    double f   = args(1).double_value();
    double x0  = args(2).double_value();
    double y0  = args(3).double_value();
    Field * field = mapToField( args(4).map_value() );
    field->t_lens(R, f, x0, y0);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

