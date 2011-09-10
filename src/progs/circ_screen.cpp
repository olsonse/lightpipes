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


#include <lightpipes/Field.h>

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {

    double r = 0.0, x0 = 0.0, y0 = 0.0;

    /*
     * Processing the command line argument 
     */

    if ( argc < 2 || argc > 4 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &r );
    if ( argc == 3 )
        sscanf ( argv[2], "%le", &y0 );
    if ( argc == 4 ) {
        sscanf ( argv[2], "%le", &y0 );
        sscanf ( argv[3], "%le", &x0 );
    }

    lightpipes::Field * field = lightpipes::Field::read (  );
    field->circular_screen( r, x0, y0 );
    field->write (  );

    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr )
{
    fprintf ( stderr, "\n%s filters the field through \
nontransparent circular screen of radius R\n", arr );

    fprintf ( stderr, "\nUSAGE: %s R [DX DY], where R \
is the radius \nof the screen in [units you use],\n\
Dx and DY are zhe shifts of the screen center\n\n", arr );


}
