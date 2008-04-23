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
    double z = 0.0, f = 0.0;

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
    sscanf ( argv[1], "%le", &f );
    sscanf ( argv[2], "%le", &z );

    if ( f == 0 )
        fprintf ( stderr, "lens_forvard: Lens with ZERO focal "
                  "length does not exist.\n" );

    Field *field = Field::read (  );
    field->lens_fresnel ( f, z );

    field->write (  );
    return EXIT_SUCCESS;
}


void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s propagates the field to "
              "distance Z [units you use]\nin VARIABLE COORDINATE  system",
              arr );
    fprintf ( stderr, "\n\nUSAGE:  " );
    fprintf ( stderr, "%s f Z, where  f is "
              "the focal lenght of the input lens,\n"
              "Z is the distance to propagate\n\n", arr );
}
