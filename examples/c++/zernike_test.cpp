#include <lightpipes/Field.h>

#include <fstream>

int main() {
  std::ofstream fint("intensity.pgm"), fpha("phase.pgm");
  using namespace lightpipes;
  Field f = Field( 256, 0.01, 1e-6 );
  f.fill( std::complex<double>(1.0,0.0) )
   .gaussian_aperture( 0.00225 )
   .zernike( 2, 2, 0.0045, 20)
   .zernike( 2, 0, 0.0045, -10)
   .fresnel(1.55)
   //.interpolate( 0, 0, 0, 0, M_PI/32 ) /* example of rotating the grid */

   /* print the field intensity and phase */
   /*            stream   out_grid   gamma max  asii? */
   .print_norm ( fint, f.info.number, 2.0, 100, true );
  f.print_phase( fpha, f.info.number, 2.0, 100, true );
  return 0;
}
