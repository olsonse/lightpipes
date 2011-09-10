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

    double f, x0, y0;

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
    sscanf ( argv[1], "%le", &f );
    x0 = y0 = 0.;
    if ( argc == 3 )
        sscanf ( argv[2], "%le", &y0 );
    if ( argc == 4 ) {
        sscanf ( argv[2], "%le", &y0 );
        sscanf ( argv[3], "%le", &x0 );
    }

    lightpipes::Field * field = lightpipes::Field::read (  );

    if ( field->info.sph_coords_factor != 0. ) {
        fprintf ( stderr, "Lens can not be applied in spherical "
                          "coordinates,\nuse CONVERT first\n" );
        exit ( 1 );
    }

    field->lens ( f, x0, y0 );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}



void error_print ( const char *arr ) {
    fprintf ( stderr,
              "\n%s filters the beam through the quadratic phase corrector\n",
              arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr, "%s  F[DX DY], where F is the focal length \n\
in [units you use],\n\
DX and DY are the transversal shifts of the lens optical axis \n\n", arr );

}
