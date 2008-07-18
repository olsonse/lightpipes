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

void error_print(const char *arr = NULL);

int main(int argc, const char *argv[]) {

    double B, dx = 0.0, dy = 0.0, angle = 0.0, mag = 1.0;
    int n_grid = 0;

    /* Processing the command line argument  */

    if (argc < 2 || argc > 7){
        error_print ( argv[0] );
        exit( 1 );
    }

    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &B );

    if ( argc >= 3 )
        sscanf ( argv[2], "%d", &n_grid );

    if ( argc >= 4 )
        sscanf ( argv[3], "%le", &dx );

    if ( argc >= 5 )
        sscanf ( argv[4], "%le", &dy );

    if ( argc >= 6 )
        sscanf ( argv[5], "%le", &angle );

    if ( argc >= 7 )
        sscanf ( argv[6], "%le", &mag );


    /* change to radians. */
    angle *= M_PI / 180.0;

    Field * field = Field::read (  );
    field->interpolate(B,n_grid,dx,dy,angle,mag);
    field->write (  );
    delete field;

    return 0;
}


void error_print ( const char *arr ) {
    fprintf ( stderr,
       "\n%s  transfers the field into a grid with different size and dimension\n",
       arr );

    fprintf ( stderr,
        "\nUSAGE: %s B [N  [X_shift  [Y_shift  [A  [M]]]]]  where \n"
        "B is the new side length of the grid; the area of the grid will"
        "    be BxB [units you use]\n"
        "    Use '-1' to indicate that the input grid size should be used\n"
        "N is the new dimension of the grid\n"
        "    Use '-1' to indicate that the input grid length should be used\n"
        "    [default:  -1 -- dimensions of input grid]\n"
        "X_shift and Y_shift are the shifts of the field in a new grid\n"
        "    [default:  0.0]\n"
        "A is the angle of rotation (after shifts) and\n"
        "    [default:  0.0]\n"
        "M (M><0) is the magnification (the last applied)\n"
        "    [default:  0.0]\n"
        "\n", arr );
}

