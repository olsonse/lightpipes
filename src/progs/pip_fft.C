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
    int ind;

    /*
     * Processing the command line argument 
     */
    if ( argc != 2 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%d", &ind );

    Field * field = Field::read (  );
    field->pip_fft ( ind );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s  performs FFT of the field distribution\n",
              arr );
    fprintf ( stderr, "\n\nUSAGE:  " );
    fprintf ( stderr, "%s IND , where IND is direction,\n"
                      "IND=1 for forward, -1 for inverse \n\n",
              arr );

}
