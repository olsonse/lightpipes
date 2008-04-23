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
#include <fstream>
#include <lightpipes/Field.h>
#define GAMMA 2.0

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {

    int output_size = 0, max_val = 255;
    bool ascii = false;

    double gamma = GAMMA;
    /*
     * Processing the command line argument 
     */
    if ( argc < 2 || argc > 6 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    if ( argc >= 3 ) {
        if ( ( strstr ( argv[2], "asc" ) ) != NULL ) {
            ascii = true;
        }
    }

    if ( argc >= 4 ) {

        if ( ( strstr ( argv[3], "sam" ) ) != NULL ) {
            output_size = 0;
        } else {
            if ( ( sscanf ( argv[3], "%d", &output_size ) ) != 0 ) {
            } else output_size = 0;
        }
    }

    if ( argc >= 5 ) sscanf ( argv[4], "%le", &gamma );
    if ( argc >= 6 ) sscanf ( argv[5], "%d", &max_val );


    std::ofstream output (argv[1]);

    Field * field = Field::read (  );
    field->print_field (output, output_size, gamma, max_val, ascii);
    output.close (  );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}




void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s writes intensity distribution\
into *.pgm file F\n", arr );

    fprintf ( stderr, "\nUSAGE: %s F ['raw'|'ascii', N, G, MAX], where F is the output filename,\n\
optional parameter N is the grid size, N=128 if omitted in command line\n\
N equals to grid sampling if you pass the word \"same\"\n\
G is the Gamma parameter, [0.1...10], higher G gives better\n\
contrast in low intensities, default G=2.0\n\
MAX is the number of gray levels, default MAX=255\n\n",
              arr );
    fprintf ( stderr,
              "Output file F can be processed with netpbm package \n\n" );


}

#undef GAMMA
