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

DEFUN_DLD (LPPipFFT, args, nargout,
"LPPipFFT\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPPipFFT performs a Fourier transform to the Field\n"
"\n"
"	USAGE:\n"
"\n"
"	Fout = LPPipFFT(Dir,Fin);\n"
"	with:\n"
"	Dir  = 1:  Forward transform\n"
"	Dir  = -1: Inverse transform\n"
"	Fin  = input field struct (see LPBegin)\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() != 2 ||
        !args(0).is_real_scalar() ||
        !args(1).is_map() || !mapIsValidField(octave_stdout, args(1).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;

        print_usage("LPPipFFT");
        return octave_value();
    }

    int dir   = args(0).int_value();
    Field * field = mapToField( args(1).map_value() );
    field->pip_fft(dir);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

