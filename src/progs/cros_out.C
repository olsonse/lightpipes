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
#include <fstream>

void error_print ( char *arr = NULL );

int main ( int argc, char *argv[] ) {
    double *int1 = NULL;
    double int2;
    double *phase1 = NULL;
    double phase2;

    int i,
        j,
        jj;
    double dx;
    long ik1;

    /*
     * Processing the command line argument 
     */

    if ( argc < 2 || argc > 2 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    Field * f = Field::read ( );
    Field & field = *f;

    dx = field.info.side_length / ( field.info.number - 1. );

    /*
     * allocating memory for output arrays 
     */

    int1   = new double[field.info.number];
    phase1 = new double[field.info.number];
    if ( int1 == NULL || phase1 == NULL) {
        fprintf ( stderr, "Allocation error in cros_out , int1\n" );
        exit ( 1 );
    }

    /*
     * writing the intensity into arrays 
     */

    i = field.info.number / 2 + 1;
    for ( j = 1; j <= field.info.number; j += 1 ) {
        ik1 = ( i - 1 ) * field.info.number + j - 1;
        jj = j - 1;
        int1[jj]   = norm(field[ik1]);
        phase1[jj] = arg(field[ik1]);

    }
#ifdef _DJGPP_
    setmode ( fileno ( stdout ), O_BINARY );

    std::ofstream fr( argv[1], "wb" );
#else
    std::ofstream fr( argv[1] );
#endif


    j = field.info.number / 2 + 1;
    for ( i = 1; i <= field.info.number; i += 1 ) {
        double cc;
        ik1 = ( i - 1 ) * field.info.number + j - 1;
        jj = i - 1;
        cc = dx * ( i - field.info.number / 2 - 1 );

        int2   = norm(field[ik1]);
        phase2 = arg(field[ik1]);

        fr << cc << '\t'
           << int1[jj] << '\t'
           << int2 << '\t'
           << phase1[jj] << '\t'
           << phase2 << '\n';
    }


    fr.close();


    field.write ( );

    return 0;
}




void error_print ( char *arr ) {
    fprintf ( stderr, "\n%s  writes X and Y cross sections of intensity \
and phase distributions into file F\n", arr );

    fprintf ( stderr, "USAGE: %s  F, where F is the output filename,\n\n\
first column   :    coordinate\n\n\
second column  :    intensity distribution along X coordinate,\n\
section through the grid center \n\
\nthird column   :    intensity distribution along Y coordinate,\n\
section through the grid center\n\
\nfourth column  :    phase distribution along X coordinate, \n\
section through the grid centern\n\
\nfifth  column  :    phase distribution along X coordinate, \n\
section through the grid center\n\
    \n", arr );
    fprintf ( stderr,
              "The distributions  can be plotted with gnuplot plot command\n\n" );


}
