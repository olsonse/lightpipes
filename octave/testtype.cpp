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
#include <octave/ov-struct.h>
#include <octave/pager.h>

DEFUN_DLD (testtype, args, nargout,
"testtype\n"
" There is no help here!\n"
" usage : testtype(map, int)\n"
"       where map is a map with a field called 'a'\n"
"*****************************************************\n"
) {
    if (args.length() != 3 ||
        !args(0).is_map() ||
        !args(1).is_real_scalar() ||
        !args(2).is_real_scalar()) {
        print_usage("testtype");
        return octave_value();
    }

    Octave_map field = args(0).map_value();

    octave_stdout << "arg(1) = " << args(1).int_value() << "\n"
                     "arg(2) = " << args(2).double_value() << std::endl;

    if (!field.contains("a")) {
        print_usage("testtype");
        return octave_value();
    }

    int val = field.contents("a")(0).int_value();
    octave_stdout << "field.a = " << val << std::endl;

    Octave_map map;
    map.assign("bob", 0);
    octave_struct retval(map);
    retval.print(octave_stdout);

    return octave_value(field);
}

