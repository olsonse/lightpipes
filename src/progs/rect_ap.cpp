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
#include <math.h>

void error_print ( const char *arr = NULL);

int main ( int argc, char *argv[] ) {

    double Lx, Ly, x0, y0, angle;

    /*
     * Processing the command line argument 
     */

    if ( argc < 2 || argc > 6 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &Ly );
    Lx = Ly;
    x0 = y0 = angle = 0.;
    if ( argc == 3 )
        sscanf ( argv[2], "%le", &Lx );
    if ( argc == 4 ) {
        sscanf ( argv[2], "%le", &Lx );
        sscanf ( argv[3], "%le", &y0 );
    }
    if ( argc == 5 ) {
        sscanf ( argv[2], "%le", &Lx );
        sscanf ( argv[3], "%le", &y0 );
        sscanf ( argv[4], "%le", &x0 );
    }
    if ( argc == 6 ) {
        sscanf ( argv[2], "%le", &Lx );
        sscanf ( argv[3], "%le", &y0 );
        sscanf ( argv[4], "%le", &x0 );
        sscanf ( argv[5], "%le", &angle );
        angle *= -M_PI / 180.0;
    }

    lightpipes::Field * field = lightpipes::Field::read (  );
    field->rectangular_aperture( Lx, Ly, x0, y0, angle );
    field->write (  );
    delete field;

    return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s filters the field \
through a square aperture \n", arr );

    fprintf ( stderr, "\nUSAGE: %s X [Y DX DY Angle], where X and Y are \n\
the sides of the  aperture in [units you use],\n\
DX and DY are the shifts of the aperure center\n\
Angle is the rotation of the aperture in degrees of arc\n\
The aperture is first shifted, then roatated\n\n", arr );
}
