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

    double phi, n_real = 1., n_imag = 0., x0 = 0., y0 = 0.;

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
    sscanf ( argv[1], "%le", &phi );

    if ( argc >= 3 )
        sscanf ( argv[2], "%le", &n_real );

    if ( argc >= 4 )
        sscanf ( argv[3], "%le", &n_imag );

    if ( argc >= 5 )
        sscanf ( argv[4], "%le", &x0 );

    if ( argc >= 6 )
        sscanf ( argv[5], "%le", &y0 );

    Field * field = Field::read (  );

    if ( field->info.sph_coords_factor != 0. ) {
        fprintf ( stderr, "Lens can not be applied in spherical "
                          "coordinates,\nuse CONVERT first\n" );
        exit ( 1 );
    }

    std::complex<double> n1(n_real,n_imag);

    field->axicon ( phi, n1, x0, y0 );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}



void error_print ( const char *arr ) {
    fprintf ( stderr,
              "\n%s adds phase for beam traversing an axicon\n"
              "\nUSAGE:  "
              "%s phi [n_real [n_imag [DX DY]]], where F is the focal length \n"
              "in [units you use],\n"
              "n_real and n_imag are the real and imaginary components of the\n"
              "index of refraction of the axicon.  [Default n = (1, 0)]\n"
              "DX and DY are the transversal shifts of the lens  optical axis \n\n",
              arr, arr );

}
