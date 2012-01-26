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

DEFUN_DLD (LPSteps, args, nargout,
"LPSteps\n"
"*****************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPSteps propagates the field a distance using the FFT method.\n"
"\n"
"	USAGE:\n"
"\n"
"	Fout = LPSteps(dz [,N [,n_filename [,k_filename [,X_filename [,dump_period]]]]],Fin);\n"
"	INPUT:\n"
"	dz          = propagation distance\n"
"	N           = Number of steps of size dz [default 1]\n"
" n_filename  = name of the file, which  contains the distribution\n"
"               of refractive index in gnuplot format\n"
"               if R equals to '' then this option is skipped\n"
" k_filename  = name of a file which  contains the two-dimensional\n"
"               distribution of the absorption coefficient\n"
"               if A equals to '' then this is skipped\n"
" X_filename  = filename where cross section of the beam\n"
"               is written at each dump_period-th step\n"
" dump_period = number of steps between dumping to X_filename\n"
"	Fin  = input field struct (see LPBegin)\n"
" OUTPUT:\n"
"	Fout = output field struct (see LPBegin)\n"
"*****************************************************\n"
) {
    if (args.length() < 2 ||
        ( args.length() > 1 && !args(0).is_real_scalar() ) ||
        ( args.length() > 2 && !args(1).is_real_scalar() ) ||
        ( args.length() > 3 && !args(2).is_string()      ) ||
        ( args.length() > 4 && !args(3).is_string()      ) ||
        ( args.length() > 5 && !args(4).is_string()      ) ||
        ( args.length() > 6 && !args(5).is_real_scalar() ) ||
        ( args.length() > 7 ) ||
        !args( args.length()-1 ).is_map() ||
        !mapIsValidField(octave_stdout, args( args.length()-1 ).map_value())
        ) {
        if ( args.length() < 2 ) octave_stdout << "requires as least two arguments!" << std::endl;
        if ( args.length() > 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() > 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() > 3 && !args(2).is_string()      ) octave_stdout << "arg 2 is invalid" << std::endl;
        if ( args.length() > 4 && !args(3).is_string()      ) octave_stdout << "arg 3 is invalid" << std::endl;
        if ( args.length() > 5 && !args(4).is_string()      ) octave_stdout << "arg 4 is invalid" << std::endl;
        if ( args.length() > 6 && !args(5).is_real_scalar() ) octave_stdout << "arg 5 is invalid" << std::endl;
        if ( args.length() > 7 ) octave_stdout << "too many arguments!" << std::endl;
        if ( !args(args.length()-1).is_map() ||
             !mapIsValidField(octave_stdout, args(args.length()-1).map_value()))
          octave_stdout << "last arg is invalid (should be Field from LPBegin)" << std::endl;

        print_usage("LPSteps");
        return octave_value();
    }

    int N = 1,
        dump_period = 1;
    std::string n_filename = "",
                k_filename = "",
                X_filename = "";
    double dz                           = args(0).double_value();
    if (args.length() > 2) N           = args(1).int_value();
    if (args.length() > 3) n_filename  = args(2).string_value();
    if (args.length() > 4) k_filename  = args(3).string_value();
    if (args.length() > 5) X_filename  = args(4).string_value();
    if (args.length() > 6) dump_period = args(5).int_value();

    Field * field = mapToField( args( args.length() -1 ).map_value() );
    field->steps(dz, N, n_filename, k_filename, X_filename, dump_period);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

