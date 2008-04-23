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
    double R, x0, y0, T;

    /*
     * Processing the command line argument 
     */
    if ( argc < 2 || argc > 5 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &R );
    x0 = y0 = 0.;
    T = 1.;
    if ( argc == 3 )
        sscanf ( argv[2], "%le", &y0 );
    if ( argc == 4 ) {
        sscanf ( argv[2], "%le", &y0 );
        sscanf ( argv[3], "%le", &x0 );
    }
    if ( argc == 5 ) {
        sscanf ( argv[2], "%le", &y0 );
        sscanf ( argv[3], "%le", &x0 );
        sscanf ( argv[4], "%le", &T );
    }

    Field * field = Field::read (  );
    field->gaussian_screen ( R, x0, y0, T );
    field->write (  );

    return 0;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s filters the beam through a gaussian screen\n",
              arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr, "%s R [DX DY T], where R is the 1/e (intensity) \n\
radius of the screen in [units you use],\n\
DX and DY are the shifts of the screen  center\n\
(1-T) gives the minimum transmission, default T=1\n\
i.e the default  minimum transmission == 0\n\n", arr );

}
