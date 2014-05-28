#include <lightpipes/Field.h>

#include <fstream>

namespace {
  static const double m=1;
  static const double nm=1e-9*m;
  static const double um=1e-6*m;
  static const double mm=1e-3*m;
  static const double cm=1e-2*m;
}

int main() {
  std::ofstream fint("intensity.pgm"), fpha("phase.pgm");
  using namespace lightpipes;
  Field f = Field( 256, 1*cm, 1*um );
  f.fill( std::complex<double>(1.0,0.0) )
   .gaussian_aperture( 2.25*mm )
   .zernike( 2, 2, 4.5*mm, 20)
   .zernike( 2, 0, 4.5*mm, -10)
   .fresnel(1.55*m)
   .interpolate( 0, 0, 0, 0, M_PI/32 ) /* example of rotating the grid */
   /* print the field intensity and phase */
   /*            stream   out_grid   gamma max  asii? */
   .print_norm ( fint, f.info.number, 2.0, 100, true );
  f.print_phase( fpha, f.info.number, 2.0, 100, true );
  return 0;
}
