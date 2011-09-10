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
#include <lightpipes/Field.h>
#include <fstream>

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {
    /*
     * Processing the command line argument 
     */
    if ( argc < 2 || argc > 2 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    lightpipes::Field * field = lightpipes::Field::read (  );
    std::ifstream field2_infile(argv[1]);
    lightpipes::Field * field2 = lightpipes::Field::read ( field2_infile );
    field2_infile.close (  );

    (*field) += (*field2);
    field->write (  );
    delete field;
    delete field2;
    return EXIT_SUCCESS;
}




void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s mixes two coherent fields:\n", arr );
    fprintf ( stderr, "\nUSAGE: %s  F < F1, where F and F1 are the input filenames\n\n", arr );
}
