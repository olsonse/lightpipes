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

#include <string>
#include <cstdio>

void error_print ( const char * arr = NULL );

int main ( int argc, char *argv[] ) {
  double step_size = 0.0;
  int N = 1, dump_period = 1;
  std::string n_filename, k_filename, X_filename;

  /*
   * Processing the command line argument 
   */
  if ( argc < 2 || argc > 7 )
    error_print( argv[0] );

  /*
   * reading the data from a command line 
   */
  sscanf ( argv[1], "%le", &step_size );
  if ( argc > 2 )
    sscanf ( argv[2], "%d", &N );

  if ( argc > 3 && std::string(argv[3]) != "void" )
    n_filename = argv[3];

  if ( argc > 4 && std::string(argv[4]) != "void" )
    k_filename = argv[4];

  if ( argc > 5 )
    X_filename = argv[5];

  if ( argc > 6 )
    sscanf ( argv[6], "%d", &dump_period );

  lightpipes::Field * field = lightpipes::Field::read (  );

  if ( field->info.sph_coords_factor != 0. ) {
    fprintf ( stderr, "%s can not be applied in spherical"
                      " coordinates,\nuse CONVERT first\n", argv[0] );
    exit( EXIT_FAILURE );
  }

  field->steps( step_size, N, n_filename, k_filename, X_filename, dump_period );
  field->write( );
  delete field;
  return EXIT_SUCCESS;
}

void error_print ( const char *arr ) {
  std::cerr << "\n"
  << arr << " propagates the field to distance N*Z [units you use]\n"
  "using finite difference method.\n\n"

  "USAGE: \n"
  << arr << " Z [N R A F M] , where:\n"
  "Z  :  step size\n"
  "N  :  the number of steps to perform [Default 1]\n"
  "R  :  name of the file, which  contains the distribution\n"
  "      of refractive index in gnuplot format\n"
  "      if R equals to < void > then this option is skipped\n"
  "A  :  name of a file which  contains the two-dimensional\n"
  "      distribution of the absorption coefficient\n"
  "      if A equals to < void > then this is skipped\n"
  "F  :  filename where cross section of the beam\n"
  "      is written at each M-th step\n\n"

  "Examples:\n"
  << arr << " 0.1 20\n"
  "    - performs 20 steps of 0.1 m in a free space, total distance of\n"
  "      propagation is  2 m \n"
  << arr << " 0.1 20 void abs_d\n"
  "    - performs 20 steps of 0.1 m in a medium absorption coefficient of\n"
  "      which is given in a file abs_d, refraction coefficient equals to\n"
  "      unity\n"
  << arr << " 0.1 20 void void file_out\n"
  "    - performs 20 steps in a free space writing cross section of the beam\n"
  "      into the file file_out, in the format:\n\n"
  "          x  I(x,0) I(0,y ) F(x,0) F(0,y) Z \n\n"
  "      where I is the intensity and F is phase-the same format as used by\n"
  "      cros_out.  Three-dimensional profile of propagating beam can be\n"
  "      plotted for example with gnuplot commands:\n\n"
  "          set par; splot 'file_out' using 1:6:2 w l;\n"
  << arr << " 0.1 200 refr_d abs_d out.dat 10\n"
  "    - performs 200 steps of 0.1 m in a medium absorption coefficient of\n"
  "      which is given in a file abs_d, refraction coefficient is given in\n"
  "      a file refr_d the profile of the beam is written into file out.dat\n"
  "      at each 10-th step\n"
  << arr << " 0.1 20 refr_d\n"
  "    - performs 20 steps of 0.1 m in a medium refraction coefficient of\n"
  "      which is given in a file refr_d without absorption\n\n"
  ;
  exit(EXIT_FAILURE);
}
