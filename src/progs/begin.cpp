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


/*
 * This filter forms the initial grid filled with unity field 
 */

#include <lightpipes/Field.h>

void error_print ( char * arr = NULL);

int main ( int argc, char *argv[] ) {

    int n_grid;
    double side_length, lambda;

    /*
     * Processing the command line argument 
     */

    if ( argc < 3 || argc > 4 ) {
        error_print ( argv[0] );
    }

    /*
     * reading the data from a command line 
     */

    sscanf ( argv[1], "%le", &side_length );
    sscanf ( argv[2], "%le", &lambda );
    n_grid = 256;
    if ( argc > 3 )
        sscanf ( argv[3], "%d", &n_grid );


    lightpipes::Field field = lightpipes::Field ( n_grid, side_length, lambda );
    field.fill( std::complex<double>(1.0,0.0) );
    field.write ( );
    /*
     * return (0); 
     */
    return 0;
}

void error_print ( char *arr ) {

    fprintf ( stderr,
              "\n%s forms the initial data structure for use "
                      "by all following programs\n"
              "\nUSAGE: begin B C [A],  where:                       \n\n"
              "A is the grid dimension, grid of AxA "
                      "points will be formed (if omitted, A=256) \n"
              "B is the side length of the grid; the area of the grid will"
                      " be BxB [units you use]\n"
              "C is the wavelength in [units you use]\n\n", arr );
    exit ( 1 );

}

