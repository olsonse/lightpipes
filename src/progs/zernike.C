/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*								*/
/*      (C) C++ port and modifications Spencer Olson 2006       */
/*      For bugs in C++ port, send bug reports to               */
/*      olsonse at umich.edu                                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include <lightpipes/Field.h>

void error_print ( const char * arr = NULL );


int main ( int argc, char *argv[] ) {
    double R, A;
    int n, m;

    /*
     * Processing the command line argument 
     */
    if ( argc != 5 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%d", &n );
    sscanf ( argv[2], "%d", &m );
    sscanf ( argv[3], "%le", &R );
    sscanf ( argv[4], "%le", &A );
    /*
     * fprintf(stderr," %d %d %e %e\n", n, m, R, A); 
     */

    Field * field = Field::read (  );
    field->zernike ( n, m, R, A );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr )
{
    fprintf ( stderr,
              "\n%s  introduces arbitrary Zernike aberration into the field distribution\n",
              arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr, "%s  n m R A, where  n and m  are the integer orders - \n\
see Born and Volf p. 465, sixth (corrected) edition, Pergamon, 1993,\n\
R is the radius at which the phase amplitude reaches A (in radians)\n\n",
              arr );

}

