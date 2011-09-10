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

DEFUN_DLD (LPGain, args, nargout,
"LPGain\n"
"****************************************************************************************************\n"
"\n"
"	LightPipes for Octave Optical Toolbox\n"
"\n"
"	LPGain simulates a laser gain sheet.\n"
"\n"
"	USAGE:\n"
"	Fout = LPGain(gain,L,Isat,Fin);\n"
"\n"
"	with:\n"
"	gain = the small signal gain\n"
"	L    = gain length\n"
"	Isat = the saturation intensity\n"
"	Fin              = the input field struct (see LPBegin)\n"
"	Fout             = the output field struct (see LPBegin)\n"
"\n"
"****************************************************************************************************\n"
) {
    if (args.length() != 4 ||
        !args(0).is_real_scalar() ||
        !args(1).is_real_scalar() ||
        !args(2).is_real_scalar() ||
        !args(3).is_map() || !mapIsValidField(octave_stdout, args(3).map_value())
        ) {
        if ( args.length() >= 1 && !args(0).is_real_scalar() ) octave_stdout << "arg 0 is invalid" << std::endl;
        if ( args.length() >= 2 && !args(1).is_real_scalar() ) octave_stdout << "arg 1 is invalid" << std::endl;
        if ( args.length() >= 3 && !args(2).is_real_scalar() ) octave_stdout << "arg 2 is invalid" << std::endl;

        print_usage("LPGain");
        return octave_value();
    }

    double gain   = args(0).double_value();
    double L      = args(1).double_value();
    double Isat   = args(2).double_value();
    Field * field = mapToField( args(3).map_value() );
    field->l_amplify(gain, L, Isat);
    octave_value retval(fieldToMap(*field));
    delete field;       field = NULL;

    return retval;
}

