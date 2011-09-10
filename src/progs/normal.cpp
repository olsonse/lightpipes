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

    if ( argc != 2 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    bool silent = true;
    if ( strstr ( argv[1], "y" ) != NULL ) silent = false;

    double norm_coeff = 0;

    lightpipes::Field * field = lightpipes::Field::read (  );
    field->normalize ( &norm_coeff );
    if (!silent) {
        std::cerr << "normal: Normalization coefficient = "
                  << norm_coeff << std::endl;
    }
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}


void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s normalizes the beam power to unity and"
                      "\noptionally prints the normalisation coefficient\n \n", arr );

    fprintf ( stderr, "\nUSAGE: %s [y,n] , if you pass y \n"
                      "then the normalization coefficient is printed to stderr \n"
                      "if you pass n - normalization is silent\n\n", arr );
}
