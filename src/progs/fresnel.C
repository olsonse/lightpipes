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
    double z;

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

    sscanf ( argv[1], "%le", &z );

    Field * field = Field::read (  );
    field->fresnel ( z );
    field->write (  );
    delete field;
    return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
    fprintf ( stderr, "\n%s propagates the field to "
                      "distance Z [units you use]\n", arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr, "%s  Z , where  Z is the distance to propagate\n"
                      "A convolution method of calculation of the \n"
                      "diffraction integral is implemented here.\n"
                      "There is no reflection an the grid borders.\n"
                      "On the same grid the function uses 8 times more RAM then forvard.\n\n",
              arr );

}

