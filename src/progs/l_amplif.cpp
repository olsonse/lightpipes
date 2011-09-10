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

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {
    double gain, i_sat, length;

    /*
     * Processing the command line argument 
     */

    if ( argc != 4 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &gain );
    sscanf ( argv[2], "%le", &length );
    sscanf ( argv[3], "%le", &i_sat );

    lightpipes::Field * field = lightpipes::Field::read (  );
    field->l_amplify ( gain, length, i_sat );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}


void error_print ( const char *arr ) {
    fprintf ( stderr,
        "\n%s filters the field through a laser medium\n", arr );

    fprintf ( stderr,
        "\nUSAGE: %s G L I_sat, where:\n"
        "         G     : gain,\n"
        "         L     : length of the medium\n"
        "         I_sat : saturation intensity of the medium\n", arr );
}

