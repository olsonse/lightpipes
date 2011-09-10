/*--------------------------------------------------------------*/
/*	(C) Gleb Vdovin 1993-1995				*/
/*	This file is a part of LightPipes beta 			*/
/*	(October 1995) distribution				*/
/*	Send bug reports to gleb@okotech.com     		*/
/*								*/
/*      (C) C++ port and modifications Spencer Olson 2006       */
/*      For bugs in C++ port, send bug reports to               */
/*      olsonse at umich.edu                                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include <math.h>
#include <lightpipes/Field.h>

void tilt (  );

void error_print ( const char *arr = NULL );

int main ( int argc, char *argv[] ) {
    double tx,
        ty;

    /*
     * Processing the command line argument 
     */
    if ( argc != 3 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */

    sscanf ( argv[1], "%le", &ty );
    sscanf ( argv[2], "%le", &tx );

    /*
     * fprintf(stderr,"%g %g ", tx,ty);
     */

    lightpipes::Field * field = lightpipes::Field::read (  );
    field->tilt ( tx, ty );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s introduces tilt into the field distribution\n",
              arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr,
              "%s X Y, where  X and Y  are the X and Y components in "
              "radians\n\n", arr );

}
