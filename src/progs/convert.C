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

void error_print ( const char *arr = NULL );

int main ( int argc, char *argv[] ) {
    /*
     * Processing the command line argument 
     */
    if ( argc != 1 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    Field * field = Field::read (  );
    field->spherical_to_normal_coords();
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s converts the field from spherical variable\n\
coordinate system into normal coordinate system", arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr, "%s\n\n", arr );
}
