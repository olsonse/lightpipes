/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*      (C) C++ port and modifications Spencer Olson 2006       */
/*      For bugs in C++ port, send bug reports to               */
/*      olsonse at umich.edu                                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include <math.h>
#include <string.h>

#include <lightpipes/Field.h>

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {
    /*
     * Processing the command line argument 
     */

    Field * field = Field::read (  );
    field->print_strehl ( std::cerr );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}


void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s: prints the general info to the stderr\n", arr );
}
