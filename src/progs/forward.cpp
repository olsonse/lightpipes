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
  namespace lp = lightpipes;
  double z = 0.;
  lp::Pair<double> new_side_length = 0.;
  lp::Pair<size_t> new_number = 0;


  lp::Field * f = lp::Field::read (  );
  lp::Field & field = *f;


  /*
   * Processing the command line argument 
   */
#ifdef _DJGPP_
  setmode ( fileno ( stdout ), O_BINARY );
#endif
  if ( argc < 2 || argc > 6 ) {
    error_print ( argv[0] );
    exit ( 1 );
  }

  /*
   * reading the data from a command line 
   */
  sscanf ( argv[1], "%le", &z );

  if ( argc >= 3 )
    sscanf ( argv[2], "%le", &new_side_length.first );
  else
    new_side_length.first = field.info.side_length.first;

  if ( argc >= 4 )
    sscanf ( argv[3], "%le", &new_side_length.second );
  else
    new_side_length.second = field.info.side_length.second;

  if ( argc >= 5 ) {
    sscanf ( argv[4], "%lu", &new_number.first );
  } else
    new_number.first = field.info.number.first;

  if ( argc >= 6 ) {
    sscanf ( argv[5], "%lu", &new_number.second );
  } else
    new_number.second = field.info.number.second;


  field.forward ( z, new_side_length, new_number );
  field.write();
  delete f;

  return EXIT_SUCCESS;
}





void error_print ( char *arr ) {
  fprintf ( stderr,
    "\n%s propagates the field to distance Z [units you use]\n", arr );
  fprintf ( stderr, "\nUSAGE:  " );
  fprintf ( stderr,
    "%s  Z [x y nx ny], where  Z is the distance to propagate\n"
    "    x,y and nx,ny are the new grid size and number of points after propagation into plane Z\n"
    "    direct calculation of the diffraction integral is implemented in forward\n\n", arr );
}
