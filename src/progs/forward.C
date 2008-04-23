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

void error_print ( char * arr = NULL);

int main ( int argc, char *argv[] ) {
    double z = 0.,
        new_side_length = 0.;
    int new_number = 0;


    Field * f = Field::read (  );
    Field & field = *f;


    /*
     * Processing the command line argument 
     */
#ifdef _DJGPP_
    setmode ( fileno ( stdout ), O_BINARY );
#endif
    if ( argc < 2 || argc > 4 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &z );

    if ( argc == 3 )
        sscanf ( argv[2], "%le", &new_side_length );
    else
        new_side_length = field.info.side_length;

    if ( argc == 4 ) {
        sscanf ( argv[2], "%le", &new_side_length );
        sscanf ( argv[3], "%d", &new_number );
    } else
        new_number = field.info.number;

    field.forward ( z, new_side_length, new_number );
    field.write();
    delete f;

    return EXIT_SUCCESS;
}





void error_print ( char *arr )
{
    fprintf ( stderr, "\n%s propagates the field to \
distance Z [units you use]\n", arr );
    fprintf ( stderr, "\nUSAGE:  " );
    fprintf ( stderr,
              "%s  Z [x1 n1], where  Z is the distance to propagate\n\
x1 and n1 are the new grid size and number of points after propagation into plane Z\n\
direct calculation of the diffraction integral is implemented in forward\n\n", arr );

}

