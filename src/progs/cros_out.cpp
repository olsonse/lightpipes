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
  namespace lp = lightpipes;

  /*
   * Processing the command line argument 
   */

  if ( argc < 2 || argc > 2 ) {
      error_print ( argv[0] );
      exit ( 1 );
  }

  lp::Field * f = lp::Field::read ( );
  lp::Field & field = *f;

  double dx = field.info.side_length.first / ( field.info.number.first - 1. );
  double dy = field.info.side_length.second / ( field.info.number.second - 1. );
  size_t i2 = field.info.number.first / 2;
  size_t j2 = field.info.number.second / 2;


#ifdef _DJGPP_
  setmode ( fileno ( stdout ), O_BINARY );

  std::ofstream fr( argv[1], "wb" );
#else
  std::ofstream fr( argv[1] );
#endif

  /*
   * writing the intensity of x-direction
   */
  for ( size_t i = 0; i < field.info.number.first; ++i ) {
    double x = dx * ( i - (double)i2 );
    const lp::Field::Pixel & px = field(i,j2);
    fr << x << '\t' << norm(px) << '\t' << arg(px) << '\n';
  }

  fr << '\n'; // blank line between x and y directions

  /*
   * writing the intensity of y-direction
   */
  for ( size_t j = 0; j < field.info.number.second; ++j ) {
    double y = dy * ( j - (double)j2 );
    const lp::Field::Pixel & px = field(i2,j);
    fr << y << '\t' << norm(px) << '\t' << arg(px) << '\n';
  }

  fr.close();

  field.write ( );
  return 0;
}




void error_print ( char *arr ) {
  fprintf(stderr,
    "\n%s  writes X and Y cross sections of intensity "
    "and phase distributions into file F\n", arr );

  fprintf(stderr,
    "USAGE: %s  F, where F is the output filename,\n\n"
    "  first column   :    coordinate\n"
    "\n"
    "  second column  :    intensity distribution along X coordinate,\n"
    "                      section through the grid center\n"
    "\n"
    "  third column   :    intensity distribution along Y coordinate,\n"
    "                      section through the grid center\n"
    "\n"
    "  fourth column  :    phase distribution along X coordinate, \n"
    "                      section through the grid centern\n"
    "\n"
    "  fifth  column  :    phase distribution along X coordinate, \n"
    "                      section through the grid center\n"
    "  \n", arr );
  fprintf(stderr,
    "The distributions  can be plotted with gnuplot plot command\n\n" );
}
