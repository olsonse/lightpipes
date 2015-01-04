

#include <lightpipes/Field.h>

#include <algorithm>
#include <fstream>
#include <vector>

#include <stdint.h>
#include <cmath>


/* FIXME:
  There appears to be quite a bit of inconsistency with the definition of the
  grid.  In some places, it appears that the grid is defined with N(xM) edges
  and in other places, it appers that the grid is defined with N(xM) nodes.  The
  two cases are most apparent in the calculation of dx(,dy) where
    - dx = side_length / number
      is when the number of edges is specified by N(xM) and
    - dx = side_length / (number - 1)
      is when the number of nodes is specified by N(xM).

  This needs to be re-looked at and perhaps compared with original code to
  ensure that things are correct and consistent.
*/


namespace lightpipes {


  Field::Info::Info()
    : number(0,0), side_length(0.,0.), lambda(0.),
      fft_level(0), sph_coords_factor(0.) { }


  bool Field::Info::compatible( const Info & that ) const {
    if (number      != that.number      ||
        side_length != that.side_length ||
        lambda      != that.lambda      ||
        sph_coords_factor != that.sph_coords_factor)
      return false;
    else
      return true;
  }


  std::istream & Field::Info::read(std::istream & in) {
    if ( !in.read((char*)&number, sizeof(number)) )
      throw std::runtime_error("Error while reading FIELD.number\n");

    if ( number.first / 2 != ( float ) ( number.first ) / 2. )
      throw std::runtime_error("Number of points must be even");
    if ( number.second / 2 != ( float ) ( number.second ) / 2. )
      throw std::runtime_error("Number of points must be even");


    if ( !in.read((char*)&side_length, sizeof(side_length)) )
      throw std::runtime_error("Error while reading FIELD.side_length");

    if ( !in.read((char*)&lambda, sizeof(lambda)) )
      throw std::runtime_error("Error while reading FIELD.lambda");

    if ( !in.read((char*)&fft_level, sizeof(fft_level)) )
      throw std::runtime_error("Error while reading FIELD.fft_level");

    if ( !in.read((char*)&sph_coords_factor, sizeof(sph_coords_factor)) )
      throw std::runtime_error("Error while reading FIELD.sph_coords_factor");

    return in;
  }


  std::ostream & Field::Info::write(std::ostream & out) {
    if ( !out.write((char*)&number, sizeof(number)) )
      throw std::runtime_error("Error while writing FIELD.number");

    if ( !out.write((char*)&side_length, sizeof(side_length)) )
      throw std::runtime_error("Error while writing FIELD.side_length");

    if ( !out.write((char*)&lambda, sizeof(lambda)) )
      throw std::runtime_error("Error while writing FIELD.lambda");

    if ( !out.write((char*)&fft_level, sizeof(fft_level)) )
      throw std::runtime_error("Error while writing FIELD.fft_level");

    if ( !out.write((char*)&sph_coords_factor, sizeof(sph_coords_factor)) )
      throw std::runtime_error("Error while writing FIELD.sph_coords_factor");

    return out;
  }







  namespace {
    const double M_2PI = 2.0 * M_PI;
    const std::complex<double> I(0.0,1.0);

    inline double fractionOf( const double & arg ) {
      return arg - std::floor(arg);
    }
  }


  int fresnl ( double xxa, double * ssa, double * cca );
  double polevl ( double x, double coef[], int N );
  double p1evl ( double x, double coef[], int N );


  Field::Field ( const Pair<size_t> & number,
                 const Pair<double> & side_length,
                 double lambda,
                 std::complex<double> init_fill,
                 int fft_level,
                 double sph_coords_factor ) : val(NULL) {
    info.number = number;
    info.side_length = side_length;
    info.lambda = lambda;
    info.fft_level = fft_level;
    info.sph_coords_factor = sph_coords_factor;
    init( init_fill );
  }


  Field::Field ( const Info & that_info,
                 std::complex<double> init_fill ) : val(NULL) {
    info = that_info;
    init( init_fill );
  }


  Field::Field (const Field & that) : val(NULL) {
    (*this) = that;
  }


  Field::~Field () {
    cleanup();
  }


  /** program read_field reads the FIELD from the standard input. */
  Field * Field::read ( std::istream & in ) throw (std::runtime_error) {
    Field::Info info;
    Field * retval = NULL;

    #ifdef _DJGPP_
      setmode ( fileno ( stdin ), O_BINARY );
    #endif

    info.read(in);
    retval = new Field (info);

    if ( !in.read((char*)retval->val, info.mem_size()) )
      throw std::runtime_error("Error while reading FIELD.val");

    return retval;
  } /* Field::read */


  /** program write_field writes field to the standard output. */
  std::ostream & Field::write(std::ostream & out) {
    #ifdef _DJGPP_
      setmode ( fileno ( stdout ), O_BINARY );
    #endif

    info.write(out);

    if ( !out.write((char*)val, info.mem_size()) )
      throw std::runtime_error("Error while writing FIELD.val\n");

    return out;
  }/* Field::write */


  void Field::base_init() {
    val = new Pixel[info.size()];
    if (val == NULL)
      throw std::runtime_error("field allocation error");

    fftw_configs[0] = fftw_plan_dft_2d(
      info.number.first, info.number.second,
      reinterpret_cast<fftw_complex*>(val),
      reinterpret_cast<fftw_complex*>(val), // operate in place
      FFTW_FORWARD, FFTW_ESTIMATE
    );

    fftw_configs[1] = fftw_plan_dft_2d(
      info.number.first, info.number.second,
      reinterpret_cast<fftw_complex*>(val),
      reinterpret_cast<fftw_complex*>(val), // operate in place
      FFTW_BACKWARD, FFTW_ESTIMATE
    );
  }

  Field & Field::operator= (const Field & that) {
    cleanup();
    this->info = that.info;
    this->base_init();
    std::copy( that.val, that.val+info.size(), val );
    return *this;
  }


  Field & Field::operator+= ( const Field & that ) {
    if (!compatible(that))
      throw std::runtime_error("cannot add fields that are not compatible");

    for ( long i = info.size()-1; i >= 0; --i ) {
      val[i] += that.val[i];
    }

    return *this;
  }



  Field & Field::circular_aperture(const double & r, const double & x0, const double & y0) {
    double r2 = SQR(r);
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    /*
     * Cutting the aperture 
     */

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x2  = SQR((i - i2) * dx - x0);
      for ( size_t j = 0; j < info.number.second; ++j ) {
        if ( (x2 + SQR((j - j2) * dy - y0)) > r2 )
          idx(i,j) = 0.0;
      }
    }

    return *this;
  }

  /*
   * f>0 for positive lens !!!!!  (Unlike the original lightpipes code)
   */
  Field & Field::lens ( const double & f, const double & x0, const double & y0 ) {
    double K = M_2PI / info.lambda;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double y = (j - j2) * dy - y0;
        double fi = -0.5 * K * ( x2 + SQR(y) ) / f;
        idx(i,j) *= exp(I * fi);
      }
    }

    return *this;
  }/* Field::lens */

  /*
   * f>0 for positive lens !!!!!  (Unlike the original lightpipes code)
   */
  Field & Field::t_lens ( const double & R, const double & f, const double & x0, const double & y0 ) {
    double K = M_2PI / info.lambda;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double y = (j - j2) * dy - y0;
        double fi = -0.5 * K * SQR( R - sqrt ( x2 + SQR(y) ) ) / f;
        idx(i,j) *= exp(I * fi);
      }
    }

    return *this;
  }


  Field & Field::lens_forvard(double f, double z) {
    const size_t sz = info.size();

    {
      double f1 = 0.;

      if ( info.sph_coords_factor != 0. )
        f1 = 1. / info.sph_coords_factor;
      else
        f1 = 10000000. * info.phy_size() / info.lambda;
      if ( ( f + f1 ) != 0. )
        f = ( f * f1 ) / ( f + f1 );
      else
        f = 10000000. * info.phy_size() / info.lambda;
    }

    double z1 = -z * f / ( z - f );

    this->forvard ( z1 );

    double ampl_scale = ( f - z ) / f;
    info.side_length.first  *= ampl_scale;
    info.side_length.second *= ampl_scale;
    info.sph_coords_factor = -1. / ( z - f );

    if ( z1 >= 0. ) {
      for ( long i = sz - 1; i >= 0; --i )
        val[i] /= ampl_scale;
    } else {
      Pixel * f1 = new Pixel[sz];
      if ( f1 == NULL ) {
        throw std::runtime_error("Allocation error, f1 in lens_forvard");
      }

      std::fill(f1, f1+sz, Pixel(0.0));

      for ( size_t i = 0; i < info.number.first; ++i ) {
        size_t i_i = info.number.first - i + 1;
        for ( size_t j = 0; j < info.number.second; ++j ) {
          size_t j_j = info.number.second - j + 1;
          f1[i*info.number.second + j] = idx(i_i,j_j) / ampl_scale;
        }
      }

      std::copy( f1, f1 + sz, val );
      delete[] f1;
    }

     return *this;
  }

  Field & Field::forvard ( const double & z ) {
    using std::abs; using std::exp;
    this->negate_alternate_elems();

    if ( z >= 0. )  this->fft3( 1 );
    else            this->fft3( -1 );

    /*
     * Spatial filter, (c) Gleb Vdovin 1986: 
     */
    { double z_abs = abs( z ) * info.lambda / 2.;
      double i2 = std::floor(info.number.first / 2);
      double j2 = std::floor(info.number.second / 2);
      double sign = (z >= 0.) ? 1. : -1.;

      for ( size_t i = 0; i < info.number.first; ++i ) {
        for ( size_t j = 0; j < info.number.second; ++j ) {
          double sw = SQR( (i - i2) / info.side_length.first )
                    + SQR( (j - j2) / info.side_length.second );
          // seems that z_abs * sw is the number of (-)2*pi phase wrappings...?
          idx(i,j) *= exp(sign * I * M_2PI * fractionOf( -z_abs * sw ));
        }
      }
    }


    if ( z >= 0. )  this->fft3( -1 );
    else            this->fft3( 1 );

    this->negate_alternate_elems();

    return *this;
  }

  Field & Field::fft3 ( int direction ) {
    fftw_execute( fftw_configs[direction > 0 ? 0 : 1] );
    (*this) *= 1./std::sqrt(info.size()); // normalize to sqrt(N*M)
    return *this;
  }

  Field & Field::circular_screen( const double & r, const double & x0, const double & y0 ) {
    double r2 = SQR(r);

    double dx = info.side_length.first / ( info.number.first );
    double dy = info.side_length.second / ( info.number.second );
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    /*
     * Cutting the aperture 
     */

    for ( size_t i = 0; i < info.number.first; ++i ) {
        double x = (i - i2) * dx - x0;
        for ( size_t j = 0; j < info.number.second; ++j ) {
            double y = (j - j2) * dy - y0;
            if ( SQR(x) + SQR(y) <= r2 ) {
                idx(i,j) = 0.0;
            }
        }
    }

    return *this;
  }

  Field & Field::rectangular_aperture ( const double & Lx,
                                              double   Ly,
                                        const double & x0,
                                        const double & y0,
                                        const double & angle ) { 
    if (Ly < 0.0) Ly = Lx;

    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    double cc = cos ( angle );
    double ss = sin ( angle );

    /*
     * Cutting the aperture 
     */

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double xx = (i - i2) * dx - x0;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double yy = (j - j2) * dy - y0;

        double x =  xx * cc + yy * ss;
        double y = -xx * ss + yy * cc;
        if ( std::abs(x) > Lx/2. || std::abs(y) > Ly/2. ) {
          idx(i,j) = 0.0;
        }
      }
    }

    return *this;
  }

  Field & Field::rectangular_screen ( const double & Lx,
                                            double   Ly,
                                      const double & x0,
                                      const double & y0,
                                      const double & angle ) { 
    if (Ly < 0.0) Ly = Lx;

    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    double cc = cos ( angle );
    double ss = sin ( angle );

    /*
     * Cutting the aperture 
     */

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double xx = (i - i2) * dx - x0;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double yy = (j - j2) * dy - y0;

        double x =  xx * cc + yy * ss;
        double y = -xx * ss + yy * cc;
        if ( std::abs(x) <= Lx/2. && std::abs(y) <= Ly/2. ) {
          idx(i,j) = 0.0;
        }
      }
    }

    return *this;
  }

  Field & Field::supergaussian_aperture ( const double & w,
                                          const int    & n,
                                          const double & x0,
                                          const double & y0,
                                          const double & A ) {
    using std::sqrt; using std::pow; using std::abs; using std::exp;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    double w2_inv = 1./ (SQR(w) * 2.);
    double Asqrt = sqrt( abs(A) );

    for ( size_t i = 0; i < info.number.first; ++i ) {
        double x2 = SQR( (i - i2) * dx - x0 );
        for ( size_t j = 0; j < info.number.second; ++j ) {
            double y = (j - j2) * dy - y0;
            double parg = pow(( x2 + SQR(y) ) * w2_inv, n);
            double cc = Asqrt * exp( - parg );
            idx(i,j) *= cc;
        }
    }
    return *this;
  }

  /*
   * 1-A is the transmission in the minimum point if the screen is modeled 
   * as a mirror then A is the maximum reflection and gauss_scr returns the
   * transmitted field 
   */
  Field & Field::supergaussian_screen ( const double & w,
                                        const int    & n,
                                        const double & x0,
                                        const double & y0,
                                        const double & A ) {
    using std::sqrt; using std::pow; using std::abs; using std::exp;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    double w2_inv = 1./SQR(w);

    for ( size_t i = 0; i < info.number.first; ++i ) {
        double x2 = SQR( (i - i2) * dx - x0 );
        for ( size_t j = 0; j < info.number.second; ++j ) {
            double y = (j - j2) * dy - y0;
            double parg = pow(( x2 + SQR(y) ) * w2_inv, n);
            double cc = sqrt( abs( 1. - A * exp ( -parg ) ) );
            idx(i,j) *= cc;
        }
    }
    return *this;
  }



  Field & Field::lens_fresnel ( const double & f_in, const double & z ) {
    double f1 = 0.;
    double f = f_in;

    if ( info.sph_coords_factor != 0. )
      f1 = 1. / info.sph_coords_factor;
    else
      f1 = 10000000. * info.phy_size() / info.lambda;


    if ( (f + f1) != 0. )
      f = (f * f1) / (f + f1);
    else
      f = 10000000. * info.phy_size() / info.lambda;

    double z1 = -z * f / (z - f);

    if ( z1 < 0. ) {
      throw std::runtime_error("lens_fresn cannot propagate behind "
                               "the focal point.  Use lens_forvard "
                               "instead." );
    }

    fresnel ( z1 );

    double ampl_scale = (f - z) / f;
    info.side_length.first *= ampl_scale;
    info.side_length.second *= ampl_scale;
    info.sph_coords_factor = -1. / ( z - f );


    double ampl_scale_inv = 1./ ampl_scale;
    for ( size_t i = 0; i < info.number.first; ++i ) {
      for ( size_t j = 0; j < info.number.second; ++j ) {
        idx(i,j) *= ampl_scale_inv;
      }
    }

    return *this;
  }

  Field & Field::spherical_to_normal_coords ( ) {
    if ( info.sph_coords_factor != 0. ) {
        this->lens ( 1. / info.sph_coords_factor, 0., 0. );
        info.sph_coords_factor = 0;
    }
    return *this;
  }


  Field & Field::forward( const double & z,
                          const Pair<double> & new_side_length,
                          const Pair<size_t> & new_number ) {
    throw std::runtime_error("not fixed for NxM yet");
#if 0
    double  dx_new, dx_old, x_new, y_new;

    int i, j, i_old, j_old;

    Field::Info newinfo;

    newinfo.side_length = new_side_length;
    newinfo.number = new_number;
    newinfo.lambda = info.lambda;
    newinfo.fft_level = info.fft_level;
    newinfo.sph_coords_factor = info.sph_coords_factor;

    Field * fieldPtr = new Field(newinfo);
    Field & field = *fieldPtr;

    /* what is R22? */
    double R22 = sqrt ( 1. / ( 2. * info.lambda * z ) );

    dx_new = newinfo.side_length    / ( newinfo.number    - 1. );
    dx_old = info.side_length / ( info.number - 1. );

    double on21 = info.number / 2 + 1;
    double nn21 = field.info.number / 2 + 1;

    /*
     * main cycle here starts 
     */

    long kk = 0,
         kk1 = 0;

    for ( i = 1; i <= field.info.number; i++ ) {
      x_new = ( i - nn21 ) * dx_new;

      for ( j = 1; j <= field.info.number; j++ ) {
        y_new = ( j - nn21 ) * dx_new;

        std::complex<double> & F = field[kk];
        F = 0.;
        //F.real() = F.imag() = 0.;
        double Fr = 0.0,
               Fi = 0.0;

        /*
         * Inner cycle here - integral itself 
         */

        kk1 = 0;
        for ( i_old = 1; i_old <= info.number; i_old++ ) {
          int io = (int)(i_old - on21);

          for ( j_old = 1; j_old <= info.number; j_old++ ) {
            int jo = (int)(j_old - on21);

            double c4c1, c2c3, c4s1,
                   s4c1, s2c3, c2s1,
                   s4c3, s2c1, c4s3,
                   s2s3, s2s1, c2s3,
                   s4s1, c4c3, s4s3, c2c1;


            double fc1(0.), fs1(0.),
                   fc2(0.), fs2(0.),
                   fc3(0.), fs3(0.),
                   fc4(0.), fs4(0.),
                   fr(0.), fi(0.);


            double P1 = R22 * ( 2 * ( dx_old * io - x_new ) + dx_old );
            double P2 = R22 * ( 2 * ( dx_old * jo - y_new ) - dx_old );
            double P3 = R22 * ( 2 * ( dx_old * io - x_new ) - dx_old );
            double P4 = R22 * ( 2 * ( dx_old * jo - y_new ) + dx_old );
            fresnl ( P1, &fs1, &fc1 );
            fresnl ( P2, &fs2, &fc2 );
            fresnl ( P3, &fs3, &fc3 );
            fresnl ( P4, &fs4, &fc4 );
            fr = 0.5 * val[kk1].real();
            fi = 0.5 * val[kk1].imag();

            c4c1 = fc4 * fc1;
            c2s3 = fc2 * fs3;
            c4s1 = fc4 * fs1;
            s4c1 = fs4 * fc1;
            s2c3 = fs2 * fc3;
            c2s1 = fc2 * fs1;
            s4c3 = fs4 * fc3;
            s2c1 = fs2 * fc1;
            c4s3 = fc4 * fs3;
            Fr       += fr * ( c2s3 + c4s1 + s4c1 + s2c3 - c2s1 - s4c3 - s2c1 - c4s3 );
          //F.real() += fr * ( c2s3 + c4s1 + s4c1 + s2c3 - c2s1 - s4c3 - s2c1 - c4s3 );

            s2s3 = fs2 * fs3;
            s2s1 = fs2 * fs1;
            c2c3 = fc2 * fc3;
            s4s1 = fs4 * fs1;
            c4c3 = fc4 * fc3;
            c4c1 = fc4 * fc1;
            s4s3 = fs4 * fs3;
            c2c1 = fc2 * fc1;
            Fr       += fi * ( -s2s3 + s2s1 + c2c3 - s4s1 - c4c3 + c4c1 + s4s3 - c2c1 );
          //F.real() += fi * ( -s2s3 + s2s1 + c2c3 - s4s1 - c4c3 + c4c1 + s4s3 - c2c1 );

            Fi       += fr * ( -c4c1 + s2s3 + c4c3 - s4s3 + c2c1 - s2s1 + s4s1 - c2c3 );
          //F.imag() += fr * ( -c4c1 + s2s3 + c4c3 - s4s3 + c2c1 - s2s1 + s4s1 - c2c3 );

            Fi       += fi * ( c2s3 + s2c3 + c4s1 + s4c1 - c4s3 - s4c3 - c2s1 - s2c1 );
          //F.imag() += fi * ( c2s3 + s2c3 + c4s1 + s4c1 - c4s3 - s4c3 - c2s1 - s2c1 );
            kk1++;
          }
        }

        F = std::complex<double>(Fr,Fi);

        kk++;
      }
    }

    (*this) = (*fieldPtr);
    delete fieldPtr;
    return *this;
#endif
  }

  Field & Field::fresnel ( const double &z ) {
    using std::sqrt;
    double dx = info.side_length.first  / info.number.first;
    double dy = info.side_length.second / info.number.second;
    if ( dx != dy )
      throw std::runtime_error("Fresnel: not fixed for non-uniform grid yet");
    size_t i2 = info.number.first / 2;
    size_t j2 = info.number.second / 2;

    /*
     * Allocating a LOT OF MEMORY 
     */

    size_t fi2 = info.number.first  * 2;
    size_t fj2 = info.number.second * 2;

    std::complex<double> * F_RI = NULL, * K_RI = NULL;

    { int len = fi2 * fj2;

      F_RI = new std::complex<double>[len];
      if ( F_RI == NULL ) {
        throw std::runtime_error("Fresnel: Allocation error, use smaller grid");
      }

      K_RI = new std::complex<double>[len];
      if ( K_RI == NULL ) {
        throw std::runtime_error("Fresnel: Allocation error, use smaller grid");
      }

      std::fill( F_RI, F_RI+len, std::complex<double>(0.0) );
      std::fill( K_RI, K_RI+len, std::complex<double>(0.0) );
    }

    /*****************************************************/

    double sq2_lz = M_SQRT2 / sqrt( info.lambda * z );
    double RRx = dx * sq2_lz;
    double RRy = dy * sq2_lz;

    int ii = 1, ij = 1;
    long ik = 0;
    const double sh = 0.5;

    for ( size_t i = info.number.first - i2; i < info.number.first + i2; ++i ) {
      double io = (double)i - info.number.first;

      double R1 = RRx * ( io - .5 + sh );
      double R3 = RRx * ( io + .5 + sh );

      for ( size_t j = info.number.second - j2; j < info.number.second + j2; ++j ) {
        double jo = (double)j - info.number.second;

        /*
         * Fresnel staff 
         * FIXME:  I suspect that the following portion including Fresnel
         * integrals is incorrect for a nonuniform grid-size.
         */
        double fc1(0.), fs1(0.),
               fc2(0.), fs2(0.),
               fc3(0.), fs3(0.),
               fc4(0.), fs4(0.);

        { double R2 = RRy * ( jo - .5 + sh );
          double R4 = RRy * ( jo + .5 + sh );

          fresnl ( R1, &fs1, &fc1 );
          fresnl ( R2, &fs2, &fc2 );
          fresnl ( R3, &fs3, &fc3 );
          fresnl ( R4, &fs4, &fc4 );
        }

        double c4c1, c2c3, c4s1,
               s4c1, s2c3, c2s1,
               s4c3, s2c1, c4s3,
               s2s3, s2s1, c2s3,
               s4s1, c4c3, s4s3, c2c1;

        c4c1 = fc4 * fc1;
        c2s3 = fc2 * fs3;
        c4s1 = fc4 * fs1;
        s4c1 = fs4 * fc1;
        s2c3 = fs2 * fc3;
        c2s1 = fc2 * fs1;
        s4c3 = fs4 * fc3;
        s2c1 = fs2 * fc1;
        c4s3 = fc4 * fs3;

        s2s3 = fs2 * fs3;
        s2s1 = fs2 * fs1;
        c2c3 = fc2 * fc3;
        s4s1 = fs4 * fs1;
        c4c3 = fc4 * fc3;
        c4c1 = fc4 * fc1;
        s4s3 = fs4 * fs3;
        c2c1 = fc2 * fc1;


        size_t ik1 = i * fj2 + j;
        double iiij = ii * ij;
        K_RI[ik1] = 0.5 * iiij * std::complex<double>(
           c4s3 + s4c3 - c4s1 - s4c1 - c2s3 - s2c3 + c2s1 + s2c1,
          -c4c3 + s4s3 + c4c1 - s4s1 + c2c3 - s2s3 - c2c1 + s2s1
        );
        /*
         * Field staff 
         */

        F_RI[ik1] = val[ik] * iiij;

        ik++;
        ij = -ij;
      }
      ii = -ii;
    }

    fftw_plan fft_plan;

    // FFTN REPLACEMENT 1
    fft_plan = fftw_plan_dft_2d(
      fi2, fj2,
      reinterpret_cast<fftw_complex*>(K_RI),
      reinterpret_cast<fftw_complex*>(K_RI), // operate in place
      FFTW_FORWARD, FFTW_ESTIMATE
    );
    fftw_execute( fft_plan );
    for ( long i = fi2*fj2 - 1; i >= 0; --i ) K_RI[i] *= 1./sqrt(fi2*fj2); // normalize
    fftw_destroy_plan(fft_plan);

    // FFTN REPLACEMENT 2
    fft_plan = fftw_plan_dft_2d(
      fi2, fj2,
      reinterpret_cast<fftw_complex*>(F_RI),
      reinterpret_cast<fftw_complex*>(F_RI), // operate in place
      FFTW_FORWARD, FFTW_ESTIMATE
    );
    fftw_execute( fft_plan );
    for ( long i = fi2*fj2 - 1; i >= 0; --i ) F_RI[i] *= 1./sqrt(fi2*fj2); // normalize
    fftw_destroy_plan(fft_plan);

    { size_t ik = 0;
      int ii = 1, ij = 1;
      for ( size_t i = 0; i < fi2; ++i ) {
        for ( size_t j = 0; j < fj2; ++j ) {
          F_RI[ik] = ((double)ii * ij) * K_RI[ik] * F_RI[ik];
          ik++;
          ij = -ij;
        }
        ii = -ii;
      }
    }

    delete[] K_RI;

    // FFTN REPLACEMENT 3
    fft_plan = fftw_plan_dft_2d(
      fi2, fj2,
      reinterpret_cast<fftw_complex*>(F_RI),
      reinterpret_cast<fftw_complex*>(F_RI), // operate in place
      FFTW_BACKWARD, FFTW_ESTIMATE
    );
    fftw_execute( fft_plan ); // don't normalize
    fftw_destroy_plan(fft_plan);


    { size_t ik = 0;
      int ii = 1, ij = 1;
      for ( size_t i = info.number.first - i2; i < info.number.first + i2; ++i ) {
        for ( size_t j = info.number.second - j2; j < info.number.second + j2; ++j ) {
          size_t ik1 =   i       * fj2 +   j      ;
          size_t ik2 = ( i - 1 ) * fj2 +   j      ;
          size_t ik3 = ( i - 1 ) * fj2 + ( j - 1 );
          size_t ik4 =   i       * fj2 + ( j - 1 );

          val[ik] = (0.25 * ii * ij)
                  * ( F_RI[ik1] - F_RI[ik2] + F_RI[ik3] - F_RI[ik4] );
          ik++;
          ij = -ij;
        }
        ii = -ii;
      }
    }

    delete[] F_RI;

    return *this;
  }




  /*
   * ----------Fresnel Integrals start here --------------- 
   */
  /*
   * Fresnel Integrals from CEPHES Modified for portability by G. Vdovin
   * (gleb@okotech.com ) Oct 1995 
   */

  /*
   * Cephes Math Library Release 2.1: January, 1989 Copyright 1984, 1987,
   * 1989 by Stephen L. Moshier 
   */



  /** Fresnel integrals.
   * PROTOTYPE:
   *   int fresnl( double x, double * S, double *C );
   * DESCRIPTION: Evaluates the Fresnel integrals
   *   x - | | C(x) = | cos(pi/2 t**2) dt, | | - 0
   *   x - | | S(x) = | sin(pi/2 t**2) dt. | | - 0
   * The integrals are approximated by rational functions 
   * if x is small. For large x, auxiliary functions f(x) and g(x) are
   * employed such that
   *           C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
   *           S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 ) .
   * ACCURACY: Relative error.
   *    Arithmetic function domain  # trials peak    rms
   *    IEEE       S(x)     0, 10   10000    2.0e-15 3.2e-16
   *    IEEE       C(x)     0, 10   10000    1.8e-15 3.3e-16
   *     DEC       S(x)     0, 10    6000    2.2e-16 3.9e-17
   *     DEC       C(x)     0, 10    5000    2.3e-16 3.9e-17 
   *
  
   * 
  Fresnel Integrals from CEPHES  
  Modified for portability by G. Vdovin (gleb@okotech.com     )
  Oct 1995 
   *
  
   *
  Cephes Math Library Release 2.1:  January, 1989
  Copyright 1984, 1987, 1989 by Stephen L. Moshier
  *
  
   * SYNOPSIS:
   *
   * double x, S, C;
   * int fresnl();
   *
   * fresnl( x, _&S, _&C );
   *
   *
   * DESCRIPTION:
   *
   * Evaluates the Fresnel integrals
   *
   *           x
   *           -
   *          | |
   * C(x) =   |   cos(pi/2 t**2) dt,
   *        | |
   *         -
   *          0
   *
   *           x
   *           -
   *          | |
   * S(x) =   |   sin(pi/2 t**2) dt.
   *        | |
   *         -
   *          0
   *
   *
   * The integrals are approximated by rational functions if x is small.
   * For large x, auxiliary functions f(x) and g(x) are employed
   * such that
   *
   * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
   * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 ) .
   *
   *
   *
   * ACCURACY:
   *
   *  Relative error.
   *
   * Arithmetic  function   domain     # trials      peak         rms
   *   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
   *   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
   *   DEC        S(x)      0, 10        6000       2.2e-16     3.9e-17
   *   DEC        C(x)      0, 10        5000       2.3e-16     3.9e-17
   */
  int fresnl ( double xxa, double * ssa, double * cca ) {
    double f, g, cc, ss, c, s, t, u;
    double x, x2;

    static double sn[6] = {
      -2.99181919401019853726E3,
      7.08840045257738576863E5,
      -6.29741486205862506537E7,
      2.54890880573376359104E9,
      -4.42979518059697779103E10,
      3.18016297876567817986E11,
    };
    static double sd[6] = {
      /*
       * 1.00000000000000000000E0,
       */
      2.81376268889994315696E2,
      4.55847810806532581675E4,
      5.17343888770096400730E6,
      4.19320245898111231129E8,
      2.24411795645340920940E10,
      6.07366389490084639049E11,
    };


    static double cn[6] = {
      -4.98843114573573548651E-8,
      9.50428062829859605134E-6,
      -6.45191435683965050962E-4,
      1.88843319396703850064E-2,
      -2.05525900955013891793E-1,
      9.99999999999999998822E-1,
    };
    static double cd[7] = {
      3.99982968972495980367E-12,
      9.15439215774657478799E-10,
      1.25001862479598821474E-7,
      1.22262789024179030997E-5,
      8.68029542941784300606E-4,
      4.12142090722199792936E-2,
      1.00000000000000000118E0,
    };

    static double fn[10] = {
      4.21543555043677546506E-1,
      1.43407919780758885261E-1,
      1.15220955073585758835E-2,
      3.45017939782574027900E-4,
      4.63613749287867322088E-6,
      3.05568983790257605827E-8,
      1.02304514164907233465E-10,
      1.72010743268161828879E-13,
      1.34283276233062758925E-16,
      3.76329711269987889006E-20,
    };
    static double fd[10] = {
      /*
       * 1.00000000000000000000E0,
       */
      7.51586398353378947175E-1,
      1.16888925859191382142E-1,
      6.44051526508858611005E-3,
      1.55934409164153020873E-4,
      1.84627567348930545870E-6,
      1.12699224763999035261E-8,
      3.60140029589371370404E-11,
      5.88754533621578410010E-14,
      4.52001434074129701496E-17,
      1.25443237090011264384E-20,
    };


    static double gn[11] = {
      5.04442073643383265887E-1,
      1.97102833525523411709E-1,
      1.87648584092575249293E-2,
      6.84079380915393090172E-4,
      1.15138826111884280931E-5,
      9.82852443688422223854E-8,
      4.45344415861750144738E-10,
      1.08268041139020870318E-12,
      1.37555460633261799868E-15,
      8.36354435630677421531E-19,
      1.86958710162783235106E-22,
    };
    static double gd[11] = {
      /*
       * 1.00000000000000000000E0,
       */
      1.47495759925128324529E0,
      3.37748989120019970451E-1,
      2.53603741420338795122E-2,
      8.14679107184306179049E-4,
      1.27545075667729118702E-5,
      1.04314589657571990585E-7,
      4.60680728146520428211E-10,
      1.10273215066240270757E-12,
      1.38796531259578871258E-15,
      8.39158816283118707363E-19,
      1.86958710162783236342E-22,
    };

    x = std::abs ( xxa );
    x2 = x * x;
    if ( x2 < 2.5625 ) {
      t = x2 * x2;
      ss = x * x2 * polevl ( t, sn, 5 ) / p1evl ( t, sd, 6 );
      cc = x * polevl ( t, cn, 5 ) / polevl ( t, cd, 6 );
      goto done;
    }
    if ( x > 36974.0 ) {
      cc = 0.5;
      ss = 0.5;
      goto done;
    }
    /*
     * Auxiliary functions for large argument 
     */
    x2 = x * x;
    t = M_PI * x2;
    u = 1.0 / ( t * t );
    t = 1.0 / t;
    f = 1.0 - u * polevl ( u, fn, 9 ) / p1evl ( u, fd, 10 );
    g = t * polevl ( u, gn, 10 ) / p1evl ( u, gd, 11 );

    t = M_PI_2 * x2;
    c = cos ( t );
    s = sin ( t );
    t = M_PI * x;
    cc = 0.5 + ( f * s - g * c ) / t;
    ss = 0.5 - ( f * c + g * s ) / t;
  done:
    if ( xxa < 0.0 ) {
      cc = -cc;
      ss = -ss;
    }
    *cca = cc;
    *ssa = ss;
    return ( 0 );
  }


  /*
   * Cephes Math Library Release 2.1: December, 1988 Copyright 1984, 1987,
   * 1988 by Stephen L. Moshier Direct inquiries to 30 Frost Street,
   * Cambridge, MA 02140 
   */


  /** Evaluate polynomial.
   * PROTOTYPE:
   *   double polevl( double x, * double coef[N+1], int N );
   * DESCRIPTION:
   * Evaluates polynomial of degree N:
   *   2 N y = C + C x + C x +...+ C x 0 1 2 N
   * Coefficients are stored in reverse order:
   *   coef[0] = C , ..., coef[N] = C .  N 0
   *
   * The function p1evl() assumes that coef[N] = 1.0 and 
   * is omitted from the array.  Its calling arguments are otherwise the
   * same as polevl().
   * SPEED:
   * In the interest of speed, there are no checks for out of bounds arithmetic.
   * This routine is used by most of the functions in the library.  Depending on
   * available equipment features, the user may wish to rewrite the program in
   * microcode or assembly language. 
   *
   *							polevl.c
   *							p1evl.c
   *
   *	Evaluate polynomial
   *
   *
   *
   * SYNOPSIS:
   *
   * int N;
   * double x, y, coef[N+1], polevl[];
   *
   * y = polevl( x, coef, N );
   *
   *
   *
   * DESCRIPTION:
   *
   * Evaluates polynomial of degree N:
   *
   *                     2          N
   * y  =  C  + C x + C x  +...+ C x
   *        0    1     2          N
   *
   * Coefficients are stored in reverse order:
   *
   * coef[0] = C  , ..., coef[N] = C  .
   *            N                   0
   *
   *  The function p1evl() assumes that coef[N] = 1.0 and is
   * omitted from the array.  Its calling arguments are
   * otherwise the same as polevl().
   *
   *
   * SPEED:
   *
   * In the interest of speed, there are no checks for out
   * of bounds arithmetic.  This routine is used by most of
   * the functions in the library.  Depending on available
   * equipment features, the user may wish to rewrite the
   * program in microcode or assembly language.
   *
   */
  double polevl ( double x, double coef[], int N ) {
    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;

    do
      ans = ans * x + *p++;
    while ( --i );

    return ( ans );
  }

  /*
   * N Evaluate polynomial when coefficient of x is 1.0.
   * Otherwise same as polevl. 
   * PROTOTYPE:
   *   double polevl( double x, * double coef[N+1], int N );
   * DESCRIPTION:
   * Evaluates polynomial of degree N:
   *   2 N y = C + C x + C x +...+ C x 0 1 2 N
   * Coefficients are stored in reverse order:
   *   coef[0] = C , ..., coef[N] = C .  N 0
   *
   * The function p1evl() assumes that coef[N] = 1.0 and 
   * is omitted from the array.  Its calling arguments are otherwise the
   * same as polevl().
   * SPEED:
   * In the interest of speed, there are no checks for out of bounds arithmetic.
   * This routine is used by most of the functions in the library.  Depending on
   * available equipment features, the user may wish to rewrite the program in
   * microcode or assembly language. 
   *
   *							p1evl()
   *                                          N
   * Evaluate polynomial when coefficient of x  is 1.0.
   * Otherwise same as polevl.
   */
  double p1evl ( double x, double coef[], int N ) {
    double ans;
    double *p;
    int i;

    p = coef;
    ans = x + *p++;
    i = N - 1;

    do
      ans = ans * x + *p++;
    while ( --i );

    return ( ans );
  }

  Field & Field::tilt ( double tx, double ty ) {
    using std::exp;
    double K = 2. * M_PI / info.lambda;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x = (i - i2) * dx;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double y = (j - j2) * dy;
        double fi = -( tx * x + ty * y ) * K;
        idx(i,j) *= exp(I * fi);
      }
    }

    return *this;
  }



  /***************************************************************/
  /*
   * Zernike polynomial
   * 
   * +-m R as in Born and Volf, p. 465, sixth edition, Pergamon n
   * 
   * The implementation have not been optimized for speed.
   * 
   */

  /* * Factorial function */
  double factorial ( int n ) throw (std::runtime_error) {
    double product;

    if ( n < 0 ) {
      throw std::runtime_error("factorial: argument is negative");
    }
    if ( n == 0 )
      return 1.;
    else {
      product = 1;
      while ( n >= 1 ) {
        product *= n;
        --n;
      }
      return product;
    }
  }

  double Zernike ( int n, int m, double rho, double phi ) throw (std::runtime_error) {
    if ( n < 0 ) {
      throw std::runtime_error(
        "Zernike: n must be >0; |m| must be less or equal than n\n"
        "if n is odd then m must be odd,\n"
        "if n is even then m must be even\n"
      );
    }

    int ind = 0;
    for ( int ncheck = n; ncheck >= -n; ncheck -= 2 ) {
      if ( ncheck == m )
        ind = 1;
    }
    if ( ind == 0 ) {
      throw std::runtime_error(
        "Zernike: n must be >0; |m| must be less or equal than n\n"
        "if n is odd then m must be odd,\n"
        "if n is even then m must be even\n"
      );
    }

    int mm = ( int ) std::abs ( m );
    double sum = 0;
    int int_sign = 1;
    for ( int s = 0; s <= ( int ) ( ( n - mm ) / 2 ); s++ ) {
      double product = 1;
      if ( n - 2 * s != 0 ) {
        product = pow ( rho, ( double ) ( n - 2 * s ) );
      } /* else product = 1; */

      product *= factorial( n - s ) * int_sign;
      product /= factorial( s )
               * factorial( ( ( n + mm ) / 2 ) - s )
               * factorial( ( ( n - mm ) / 2 ) - s );
      sum += product;
      int_sign = -int_sign;
    }

    if ( m <= 0 )
      return sum * cos ( m * phi );
    else
      return sum * sin ( m * phi );
  }
/*****************************************************************/
/*****************        END of Zernike     *********************/


  Field & Field::zernike ( int n, int m, double R, double A ) {
    using std::sqrt; using std::arg; using std::exp;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x = (i - i2) * dx;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double y = (j - j2) * dy;
        double rho = sqrt( SQR(x) + SQR(y) ) / R;
        double phi = arg( std::complex<double>(x, y) );

        double fi = A * Zernike ( n, m, rho, phi );
        idx(i,j) *= exp(I * fi);
      }
    }

    return *this;
  }

  std::pair<double,double> Field::get_strehl_and_energy () const {
    using std::norm; using std::sqrt;
    /*
     * Calculating the power 
     */
    double sum = 0., sum1r = 0., sum1i = 0., sum2 = 0.;
    for ( size_t i = 0; i < info.number.first; ++i ) {
      for ( size_t j = 0; j < info.number.second; ++j ) {
        const Pixel & px = idx(i,j);
        double s_p = norm(px);
        sum2  += s_p;
        sum   += sqrt( s_p );
        sum1r += px.real();
        sum1i += px.imag();
      }
    }
    double sum1 = SQR(sum1r) + SQR(sum1i);

    if ( sum == 0 )
      throw std::runtime_error("Strehl: Zero beam power, program terminated");

    return std::make_pair(sum1 / SQR(sum), sum2);
  }

  double Field::get_strehl () const {
    return this->get_strehl_and_energy().second;
  }

  std::ostream & Field::print_strehl (std::ostream & output) const {
    using std::norm;
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    { Pair<double> strehl_energy = get_strehl_and_energy();
      output << "Strehl: ratio= " << strehl_energy.first
             << " energy= " << (strehl_energy.second * (dx*dy))
             << std::endl;
    }

    /*
     * Calculating the center of gravity: 
     */
    double sum = 0., sum1r = 0., sum1i = 0.;
    for ( size_t i = 0; i < info.number.first; ++i ) {
      double y = (i - i2) * dx;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double x = (j - j2) * dy;
        double sum2 = norm( idx(i,j) );
        sum1r += sum2 * x;
        sum1i += sum2 * y;
        sum += sum2;
      }
    }

    const double x_c = sum1r / sum;
    const double y_c = sum1i / sum;

    output << "Center_of_gravity: x= " << x_c << " y= " << y_c << std::endl;

    /*
     * Calculating moments of the distribution 
     */
    double sum1x = 0., sum1y = 0.;

    sum1r = 0.;
    for ( size_t i = 0; i < info.number.first; ++i ) {
      double y = (i - i2) * dx;
      double y_y_c = y - y_c;
      for ( size_t j = 0; j < info.number.second; ++j ) {
        double x = (j - j2) * dy;
        double x_x_c = x - x_c;
        double intens = norm( idx(i,j) );
        sum1r += intens * ( SQR(x_x_c) + SQR(y_y_c) );
        sum1x += intens * SQR(x_x_c);
        sum1y += intens * SQR(y_y_c);
      }
    }

    output << "Standard deviation:  "
              " S_r= " << sqrt ( sum1r / sum )
           << " S_x= " << sqrt ( sum1x / sum )
           << " S_y= " << sqrt ( sum1y / sum ) << '\n'
           << "Grid side_length: " << info.side_length
           << ", Grid sampling: " << info.number
           << std::endl;

    return output;
  }

  Field & Field::pip_fft(const int & direction) {
    info.fft_level += direction;
    if ( info.fft_level != 0 ) this->negate_alternate_elems();

    fft3 ( direction );

    if ( info.fft_level == 0 ) this->negate_alternate_elems();

    return *this;
  }

  std::ostream & Field::print_field(std::ostream  & output,
                                    Pair<size_t>    output_size /* = 0 */,
                                    const double  & gamma /* = 2.0 */,
                                    const int     & max_val /* = 255 */,
                                    const bool    & ascii /* = false */,
                                    const bool    & output_norm /* = true. */
                                    ) {
    using std::ceil; using std::floor; using std::pow; using std::arg;
    if (output_size.first  == 0 || output_size.first > info.number.first)
      output_size.first  = info.number.first;
    if (output_size.second == 0 || output_size.second > info.number.second)
      output_size.second = info.number.second;

    size_t istep = 1, jstep = 1;

    if ( info.number.first / output_size.first > 1 ) {
      istep = info.number.first / output_size.first;
      output_size.first = (int)ceil((double)info.number.first / (double)istep);
    }

    if ( info.number.second / output_size.second > 1 ) {
      jstep = info.number.second / output_size.second;
      output_size.second= (int)ceil((double)info.number.second / (double)jstep);
    }

    /*
     * header of the PNM file 
     * seo:  
     * FIXME:  This header is wrong in size for istep!=1 and should we
     * actually use P5?
     */
    output << (ascii ? "P2\n" : "P5\n")
           << "#Creator: LightPipes (C) 1993-1996, Gleb Vdovin\n"
           << "#Creator: LightPipes/C++ (C) 2006-2008, Spencer E. Olson\n"
           << (istep==1 ? output_size.first  : (output_size.first -1)) << ' '
           << (jstep==1 ? output_size.second : (output_size.second-1)) << '\n'
           << max_val << '\n';

    /* determine the maximum intensity in the image so as to normalize to
     * max_val when writing to file. */
    double max_px = M_PI;
    if (output_norm) {
      for ( long i = info.size()-1; i >= 0; --i ) {
        double sum = norm( val[i] );
        if ( sum > max_px ) max_px = sum;
      }
    }

    /* now write pixel (set) values out to file. */
    int i_i = 1;
    for ( size_t i = 0; i <= info.number.first - istep; i += istep ) {
      for ( size_t j = 0; j <= info.number.second - jstep; j += jstep ) {
        /* determine the average intensity/phase in the next pixel set. */

        double aval = 0.0;
        if (output_norm) {
          double tot = 0;
          for ( size_t ii = i; ii < i + istep; ++ii )
            for ( size_t jj = j; jj < j + jstep; ++jj )
              tot += norm( idx(ii,jj) );

          aval = tot / (istep * jstep);
        } else {
          Pixel tot = 0;
          for ( size_t ii = i; ii < i + istep; ++ii )
            for ( size_t jj = j; jj < j + jstep; ++jj )
              tot += idx(ii,jj);

          aval = arg(tot);
        }

        /* determine the integer bitmap value. */
        int i0 = (int)floor( pow(aval/max_px, 1./(gamma + 1e-4)) * max_val );

        if (ascii) {
          output << i0 << ' ';
          ++i_i;
          if ( i_i == 40 ) {
            output << '\n';
            i_i = 1;
          }

        } else {
          if (max_val <= 255) {
            uint8_t i0T = i0;
            output.write( reinterpret_cast<char*>(&i0T), sizeof(i0T) );
          } else if (max_val <= 65535) {
            uint16_t i0T = i0;
            output.write( reinterpret_cast<char*>(&i0T), sizeof(i0T) );
          } else {
            output.write( reinterpret_cast<char*>(&i0), sizeof(i0) );
          }
        }
      } /* j loop */
    } /* i loop */
    return output;
  }


  Field & Field::normalize(double * norm_coeff /* = NULL */) {
    /*
     * Calculating the power 
     */
    double sum = 0;
    for ( long i = info.size()-1; i >= 0; --i )
      sum += norm( val[i] );
    sum *= info.side_length.first  / info.number.first
         * info.side_length.second / info.number.second;

    if ( sum == 0 )
      throw std::runtime_error("normal: Zero beam power, program terminated");

    if ( norm_coeff != NULL )
      (*norm_coeff) = sum;

    /*
     * Normalizing the power 
     */
    sum = sqrt( 1. / sum );
    for ( long i = info.size()-1; i >= 0; --i )
      val[i] *= sum;

    return *this;
  }

  Field & Field::l_amplify(const double & gain,
                           const double & length,
                           const double & i_sat) {
    using std::norm; using std::exp;
    for ( long i = info.size()-1; i >= 0; --i ) {
      Pixel & px = val[i];
      double intensity = norm( px );
      double ss = exp( length * ( gain / ( 1. + ( intensity / i_sat ) ) ) );
      px *= ss;
    }

    return *this;
  }

  Field & Field::axicon(const double & phi,
                        const std::complex<double> & n1,
                        const double & x0,
                        const double & y0 ) {
    using std::tan; using std::abs; using std::sqrt; using std::exp;
    double KtanB = (M_2PI / info.lambda) * tan( 0.5 * (M_PI - abs(phi)) );
    std::complex<double> epsilon = 1.0 - n1;
    double dx = info.side_length.first  / info.number.first;
    double dy = info.side_length.second / info.number.second;
    double i2 = std::floor(info.number.first / 2);
    double j2 = std::floor(info.number.second / 2);

    for ( size_t i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);

      for ( size_t j = 0; j < info.number.second; ++j ) {
        double r = sqrt(x2 + SQR((j - j2) * dy - y0));
        std::complex<double> fi = KtanB * ( epsilon*r );
        idx(i,j) *= exp(I * fi);
      }
    }

    return *this;
  }/* Field::axicon */



  /* *** STUFF FOR INTERPOLATE BEGINS HERE *** */
  template <class T, unsigned int N>
  class SquareMatrix {
  public:
    T val[N][N];

    void zero() {
      std::fill( val[0], val[0]+(N*N), T(0));
    }

    T & operator()(const unsigned int & i, const unsigned int & j) {
      return val[i][j];
    }

    const T & operator()(const unsigned int & i, const unsigned int & j) const {
      return val[i][j];
    }
  };




  template <class T>
  T int4 ( const double & x2, const double & dx,
           const T & y1, const T & y2, const T & y3, const T & y4,
           const double & xz );

  template <class T>
  T int16 ( double xc, double yc, double xp, double yp, double dx,
            const SquareMatrix<T,4> & z ) {
    T zz1  = int4 ( yc, dx, z(0,0), z(0,1), z(0,2), z(0,3), yp );
    T zz2  = int4 ( yc, dx, z(1,0), z(1,1), z(1,2), z(1,3), yp );
    T zz3  = int4 ( yc, dx, z(2,0), z(2,1), z(2,2), z(2,3), yp );
    T zz4  = int4 ( yc, dx, z(3,0), z(3,1), z(3,2), z(3,3), yp );
    /** FIXME:  the original code had zz1,zz2,zz2,zz3.
     * I think that was a bug.  Work it out some time and see if the following
     * is now correct. */
    return   int4 ( xc, dx, zz1,    zz2,    zz3,    zz4, xp );
  }

  /*
   * cubic four-point interpolation
   * 
   * a+b*x1+c*x1^2+d*x1^3=y1 ...................... a+b*x4+c*x4^2+d*x4^3=y4
   * 
   * where x1=x2-dx; x3=x2+dx; x4=x2+2*dx; the grid is uniform
   * 
   */

  template <class T>
  T int4 ( const double & x2, const double & dx,
           const T & y1, const T & y2, const T & y3, const T & y4,
           const double & xz ) {
    double t0;
    T a, b, c, d;

    /*
     * if(xz < x2 || xz >x2+dx ){ fprintf(stderr,"out of range in the
     * cubic interpolation routine \n"); }
     */

    /*
     * maple staff, sorry 
     */
    t0 = 1. / ( dx * dx * dx );

    a =  t0 * (  y2 * (3.0 * x2 * dx * dx) + y1 * (      x2 * x2 * x2) +
                 y1 * (3.0 * x2 * x2 * dx) + y1 * (2.0 * x2 * dx * dx) -
                 y4 * (      x2 * x2 * x2) + y3 * (3.0 * x2 * x2 * x2) +
                 y3 * (3.0 * x2 * x2 * dx) - y3 * (6.0 * x2 * dx * dx) +
                 y4 * (      x2 * dx * dx) + y2 * (6.0 * dx * dx * dx) -
                 y2 * (3.0 * x2 * x2 * x2) - y2 * (6.0 * x2 * x2 * dx) ) / 6.0;

    b = -t0 * (  y3 * (-6.0 * dx * dx) + y4 * (      dx * dx) +
                 y2 * ( 3.0 * dx * dx) + y3 * (9.0 * x2 * x2) -
                 y4 * ( 3.0 * x2 * x2) - y2 * (9.0 * x2 * x2) -
                 y2 * (12.0 * x2 * dx) + y3 * (6.0 * x2 * dx) +
                 y1 * ( 3.0 * x2 * x2) + y1 * (6.0 * x2 * dx) +
                 y1 * ( 2.0 * dx * dx) ) / 6.0;
                 
    d = -t0 * ( -3.0 * y2 - y4 + y1 + 3.0 * y3 ) / 6.0;

    c = t0 *  (  y2 * (-3.0 * x2) +
                 y3 * dx +
                 y1 * dx -
                 y2 * (2.0 * dx) -
                 y4 * x2 +
                 y3 * (3.0 * x2) +
                 y1 * x2
              ) * 0.5;

    /*
     * return a+b*xz+c*xz*xz+d*xz*xz*xz;
     */
    return a + xz * ( b + xz * ( c + xz * d ) );
  }

  template <class T>
  T inv_squares ( double x, double y, double dx,
                  const T & z,  const T & zx,
                  const T & zy, const T & zxy,
                  double x1, double y1 ) {
    using std::abs;
    double tol = 1e-6 * dx;
    if ( x1 < x - tol || x1 > x + dx + tol || y1 < y - tol
         || y1 > y + dx + tol ) {
        char err[256] = "";
        snprintf ( err, 256,
                  "out of range in inv_squares %g %g %g %g %g %g %g\n", x,
                  x1, y, y1, dx, y1 - y, x1 - x );
        throw std::runtime_error(err);
    }

    double xlow = x1 - x;
    double xhigh = x + dx - x1;
    double ylow = y1 - y;
    double yhigh = y + dx - y1;

    if ( xlow < -tol || xhigh < -tol || ylow < -tol || yhigh < -tol ) {
        char err[256] = "";
        snprintf ( err, 256,
                  " inv_squares: out of range, %g %g %g %g\n",
                  xlow, xhigh, ylow, yhigh );
        throw std::runtime_error(err);
    }

    if ( abs ( xlow  ) < tol ) return z  + ylow * ( zy - z )   / dx;
    if ( abs ( ylow  ) < tol ) return z  + xlow * ( zx - z )   / dx;
    if ( abs ( xhigh ) < tol ) return zx + ylow * ( zxy - zx ) / dx;
    if ( abs ( yhigh ) < tol ) return zy + xlow * ( zxy - zx ) / dx;

    double s1  = 1. / ( xlow * ylow );
    double s2  = 1. / ( xhigh * ylow );
    double s3  = 1. / ( xlow * yhigh );
    double s4  = 1. / ( xhigh * yhigh );
    double sum = s1 + s2 + s3 + s4;
    s1 = s1 / sum;
    s2 = s2 / sum;
    s3 = s3 / sum;
    s4 = s4 / sum;

    return z * s1 + zx * s2 + zy * s3 + zxy * s4;
  }



  Field & Field::interpolate(const Pair<double> & new_side_length /* = 0.0 */,
                             const Pair<size_t> & new_number /* = 0.0 */,
                             const double & x_shift /* = 0.0 */,
                             const double & y_shift /* = 0.0 */,
                             const double & angle /* = 0.0 */,
                             const double & magnif /* = 1.0 */) 
                             throw (std::runtime_error) {
    Info new_info = info;
    if (new_side_length.first > 0) new_info.side_length.first = new_side_length.first;
    if (new_side_length.second > 0) new_info.side_length.second = new_side_length.second;
    if (new_number.first > 0) new_info.number.first = new_number.first;
    if (new_number.second > 0) new_info.number.second = new_number.second;

    if ( magnif <= 0.0 )
      throw std::runtime_error("Field::interpolate:  magnification <= 0.0");

    Field new_field(new_info);

    SquareMatrix<std::complex<double>,4> z; z.zero();

    throw std::runtime_error("not fixed for NxM yet");
#if 0
    long n_old_max = info.number * ( info.number - 1 ) - 1;

    double dx_new = new_field.info.side_length / ( new_field.info.number - 1. );
    double dx_old = info.side_length / ( info.number - 1. );


    double on21 = info.number / 2 + 1;
    double nn2  = new_field.info.number / 2;
    double lower = ( 1 - on21 ) * dx_old;
    double upper = ( info.number - on21 ) * dx_old;

    double cc = cos ( angle );
    double ss = sin ( angle );


    for ( int i = 0; i < new_field.info.number; i++ ) {
      for ( int j = 0; j < new_field.info.number; j++ ) {

        double x0 = ( i - nn2  ) * dx_new - x_shift;
        double y0 = ( j - nn2  ) * dx_new - y_shift;
        double x_new = ( x0 * cc + y0 * ss ) / magnif;
        double y_new = ( -x0 * ss + y0 * cc ) / magnif;

        int i_old = ( int ) ( floor ( x_new / dx_old ) + on21 );
        double x_old = ( i_old - on21 ) * dx_old;
        int j_old = ( int ) ( floor ( y_new / dx_old ) + on21 );
        double y_old = ( j_old - on21 ) * dx_old;

        int i_small = 0;
        int i_local = 0;
        /*
         * first row 
         */
        long n_tmp = ( i_old - 2 ) * info.number + j_old - 2;

        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(0,0) = val[n_tmp];
        else
          z(0,0) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(1,0) = val[n_tmp];
        else
          z(1,0) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
            i_small = 1;
            i_local = 1;
        }
        if ( i_local != 1 )
          z(2,0) = val[n_tmp];
        else
          z(2,0) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(3,0) = val[n_tmp];
        else
          z(3,0) = 0.0;
        i_local = 0;

        /*
         * second row 
         */
        n_tmp = ( i_old - 2 ) * info.number + j_old - 1;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(0,1) = val[n_tmp];
        else
          z(0,1) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(1,1) = val[n_tmp];
        else
          z(1,1) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(2,1) = val[n_tmp];
        else
          z(2,1) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(3,1) = val[n_tmp];
        else
          z(3,1) = 0.0;
        i_local = 0;

        /*
         * third row 
         */
        n_tmp = ( i_old - 2 ) * info.number + j_old;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(0,2) = val[n_tmp];
        else
          z(0,2) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(1,2) = val[n_tmp];
        else
          z(1,2) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(2,2) = val[n_tmp];
        else
          z(2,2) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(3,2) = val[n_tmp];
        else
          z(3,2) = 0.0;
        i_local = 0;


        /*
         * fourth row 
         */
        n_tmp = ( i_old - 2 ) * info.number + j_old + 1;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(0,3) = val[n_tmp];
        else
          z(0,3) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(1,3) = val[n_tmp];
        else
          z(1,3) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(2,3) = val[n_tmp];
        else
          z(2,3) = 0.0;
        i_local = 0;

        n_tmp += info.number;
        if ( n_tmp < 0 || n_tmp > n_old_max ) {
          i_small = 1;
          i_local = 1;
        }
        if ( i_local != 1 )
          z(3,3) = val[n_tmp];
        else
          z(3,3) = 0.0;
        i_local = 0;




        /* now finally create the new value. */
        if ( i_small == 1 ) {
          if (    x_new > lower
               && x_new < upper
               && y_new > lower
               && y_new < upper ) {
            new_field(i,j) =
                inv_squares ( x_old, y_old, dx_old,
                              z(1,1), z(2,1),
                              z(1,2), z(2,2),
                              x_new, y_new )
                / magnif;
          }
        } else {
          if (    x_new > lower
               && x_new < upper
               && y_new > lower
               && y_new < upper ) {
            new_field(i,j) =
                int16 ( x_old, y_old, x_new, y_new, dx_old, z ) / magnif;
          }
        }
      }/* j */
    }/* i */


    /* now move the new info+data into the proper storage location. */
    *this = new_field;

    return *this;
#endif
  }

  /* *** STUFF FOR INTERPOLATE ENDS HERE *** */

  namespace {
    struct Eliminator {
      typedef std::complex<double> Pixel;

      const size_t N;
      std::vector<Pixel> alpha,
                         beta;
      const double a_i, b_i;

      Eliminator( const size_t & N, const double & dx2) :
        N(N),
        alpha(N-1, Pixel(0,0)),
        beta (N-1, Pixel(0,0)),
        a_i( -1 / dx2 ),
        b_i( -1 / dx2 )
        { }

      void operator()(const std::vector<Pixel> c,
                      const std::vector<Pixel> p,
                            std::vector<Pixel> u ) {
        /*
         * initial condition, everything is going to be zero at the edge
         */
        alpha[N-2] = beta[N-2] = alpha[0] = beta[0] = std::complex<double>(0,0);

        /*
         * forward elimination
         */
        for ( size_t i = 0; i < N - 3; ++i ) {
          std::complex<double> cc = c[i+1] -  a_i * alpha[i];
          alpha[i+1] = b_i / cc;
          beta[i+1] = ( p[i] + a_i * beta[i] ) / cc;
        }

        /* This routine seems somewhat hackish with the absorption.  Also skips
         * c[N-2] entirely. */

        /*
         * edge amplitude =0
         */
        /* Probably not important:  I think c[N-2] instead of original c[N-1] */
        u[N-1] = ( p[N-3] + a_i * beta[N-2] )
               / ( c[N-2] - a_i * alpha[N-2] );
        /*
         * backward elimination
         */
        for ( int i = N - 2; i >= 0; --i ) {
          u[i] = alpha[i] * u[i+1] + beta[i];
        }
      }
    };
  }

  /* *** STUFF FOR STEPS BEGINS HERE *** */
  Field & Field::steps( const double & step_size,
                        const int & N,
                        const std::string & n_filename,
                        const std::string & k_filename,
                        const std::string & dump_filename,
                        const int & dump_period)
                        throw (std::runtime_error) {
    Pixel * n = new Pixel[ info.size() ];
    if ( n == NULL )
      throw std::runtime_error("Field::step: Allocation error for 'n'");
    std::fill( n, n+info.size(), Pixel(1.0, 0.0) );

    if ( n_filename != "" ) {
      std::ifstream n_file( n_filename.c_str() );
      if ( !n_file )
        throw std::runtime_error(
          "Field::steps: error opening file refractive index file (%s)"
        );

      size_t ik = 0;
      for ( size_t i = 0; i < info.number.first; ++i ) {
        for ( size_t j = 0; j < info.number.second; ++j ) {
          n_file >> n[ik++].real();
          if ( !n_file )
            throw std::runtime_error(
              "Field::steps: reading the refractive indices "
              "end of input file reached, exiting"
            );
        }
      }
    }

    if ( k_filename != "" ) {
      std::ifstream k_file( k_filename.c_str() );
      if ( !k_file )
        throw std::runtime_error(
          "Field::steps: error opening file absorption index file (%s)"
        );

      size_t ik = 0;
      for ( size_t i = 0; i < info.number.first; ++i ) {
        for ( size_t j = 0; j < info.number.second; ++j ) {
          k_file >> n[ik++].imag();
          if ( !k_file )
            throw std::runtime_error(
              "Field::steps: reading the absorption indices "
              "end of input file reached, exiting"
            );
        }
      }
    }

    this->steps( step_size, N, n, dump_filename, dump_period );
    delete[] n;
    return *this;
  }

  namespace {
    typedef std::complex<double> Pixel;
    struct GetN {
      const Pixel * n;
      GetN( const Pixel * n ) : n(n) { }
      const Pixel & operator[](const size_t & i) const { return n[i]; }
    };

    static const std::complex<double> One(1.0, 0.0);
    struct VacuumN {
      const Pixel & operator[](const size_t &) const { return One; }
    };
  }/* anonymous */

  Field & Field::steps( const double & step_size,
                        const int & N,
                        const Pixel * n,
                        const std::string & dump_filename,
                        const int & dump_period)
                        throw (std::runtime_error) {
    if ( n != NULL ) {
      GetN getn(n);
      this->stepsImpl( step_size, N, getn, dump_filename, dump_period );
    } else {
      VacuumN vacuum_n;
      this->stepsImpl( step_size, N, vacuum_n, dump_filename, dump_period );
    }
    return *this;
  }

  template < typename n_accessor >
  void Field::stepsImpl(const double & step_size,
                        const int & N,
                        const n_accessor & n,
                        const std::string & dump_filename,
                        const int & dump_period)
                        throw (std::runtime_error) {
    if ( info.sph_coords_factor != 0. ) {
      throw std::runtime_error(
        "steps:  spherical coords not supported.  First call "
        "Field::spherical_to_normal_coords()"
      );
    }
    const double K = 2. * M_PI / info.lambda;
    const double dx2 = SQR( info.side_length.first / (info.number.first - 1.0 ) );
    const double dy2 = SQR( info.side_length.second / (info.number.second - 1.0 ) );

    //if ( dx2 != dy2 )
    //  throw std::runtime_error("steps: not fixed for non-uniform grid yet");

    Eliminator elim_i(info.number.first, dx2),
               elim_j(info.number.second, dy2);

    std::vector<Pixel> ui(info.number.first,    Pixel(0,0)),
                       uj(info.number.second,   Pixel(0,0)),
                       ci(info.number.first,    Pixel(0,0)),
                       cj(info.number.second,   Pixel(0,0)),
                       pi(info.number.first-2,  Pixel(0,0)),
                       pj(info.number.second-2, Pixel(0,0));

    std::ofstream dump_file;
    if ( dump_filename != "" ) {
      dump_file.open( dump_filename.c_str() );
    }

    /*
     * absorption at the borders is described here
     */
    const double AA = -10. / (step_size/2.0) / N; /* total absorption */
    const double band_pow = 2.; /* profile of the absorption border,
                                 * 2=quadratic */
    /*
     * width of the absorption border
     */
    const size_t i_left  = 0.1 * info.number.first  + 1;
    const size_t i_right = 0.9 * info.number.first  + 1;
    const size_t j_left  = 0.1 * info.number.second + 1;
    const size_t j_right = 0.9 * info.number.second + 1;


    /*
     * Main loop, steps here
     */
    for ( int istep = 1; istep <= N; ++istep ) {
      std::fill( uj.begin(), uj.end(), Pixel(0,0) );

      /*
       * Elimination in the direction i, halfstep
       */
      for ( long i = info.size() - 1; i >= 0; --i ) {
        register double fi = 0.25 * K * (step_size/2.0) * ( n[i].real() - 1. );
        val[i] *= exp(I * fi);
      }

      for ( size_t ii = 1; ii < info.number.first - 2; ii += 2 ) {
        size_t i = ii;

        for ( size_t j = 1; j < info.number.second - 1; ++j ) {
          const Pixel & uij   = idx(i    , j);
          const Pixel & ui1j  = idx(i + 1, j);
          const Pixel & ui_1j = idx(i - 1, j);

          pj[j-1] = 2. * I * K / (step_size/2.0) * uij
                  - (-2. * uij + ui1j + ui_1j) / dx2;
        }

        for ( size_t j = 0; j < info.number.second; ++j ) {
          const size_t ik = i * info.number.second + j;

          cj[j] = -2/dx2 + I * K * ( 2/(step_size/2) - n[ik].imag() );
          /*
           * absorption borders are formed here 
           */
          if ( j <= j_left )
            cj[j].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(j_left - j) /
                                     static_cast<double>(j_left),
                                     band_pow );

          if ( j >= j_right )
            cj[j].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(j - j_right) /
                                     static_cast<double>(info.number.second - j_right),
                                     band_pow );
        }

        for ( size_t j = 0; j < info.number.second; ++j ) {
          idx(i-1,j) = uj[j];
        }

        elim_j(cj, pj, uj);

        i = ii + 1;

        for ( size_t j = 1; j < info.number.second - 1; ++j ) {
          const Pixel & uij   = idx(i    , j);
          const Pixel & ui1j  = idx(i + 1, j);
          const Pixel & ui_1j = idx(i - 1, j);

          pj[j-1] = 2. * I * K / (step_size/2.0) * uij
                  - (-2. * uij + ui1j + ui_1j) / dx2;
        }

        for ( size_t j = 0; j < info.number.second; ++j ) {
          const size_t ik = i * info.number.second + j;

          cj[j] = -2/dx2 + I * K * ( 2/(step_size/2) - n[ik].imag() );
          /*
           * absorption borders are formed here 
           */
          if ( j <= j_left )
            cj[j].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(j_left - j) /
                                     static_cast<double>(j_left),
                                     band_pow );

          /* seems like this extra factor of 2 here is not consistent. */
          if ( j >= j_right )
            cj[j].imag() -= ( AA * 2. * K )
                         * std::pow( static_cast<double>(j - j_right) /
                                     static_cast<double>(info.number.second - j_right),
                                     band_pow );
        }

        for ( size_t j = 0; j < info.number.second; ++j ) {
          idx(i-1,j) = uj[j];
        }

        elim_j(cj, pj, uj);

      }

      /* something to the last row. */
      for ( size_t j = 0; j < info.number.second; ++j ) {
        idx(info.number.first-1,j) = uj[j];
      }

      for ( long i = info.size() - 1; i >= 0; --i ) {
        register double fi = 0.5 * K * (step_size/2.0) * ( n[i].real() - 1. );
        val[i] *= exp(I * fi);
      }


      /*
       * Elimination in the j direction is here, halfstep 
       */
      std::fill( ui.begin(), ui.end(), Pixel(0,0) );
      
      for ( size_t jj = 1; jj < info.number.second - 2; jj += 2 ) {
        size_t j = jj;
        for ( size_t i = 1; i < info.number.first - 1; ++i ) {
          const Pixel & uij   = idx(i,j    );
          const Pixel & uij1  = idx(i,j + 1);
          const Pixel & uij_1 = idx(i,j - 1);

          pi[i-1] = 2. * I * K / (step_size/2.0) * uij
                  - (-2. * uij + uij1 + uij_1) / dy2;
        }

        for ( size_t i = 0; i < info.number.first; ++i ) {
          const size_t ik = i * info.number.second + j;

          ci[i] = -2/dy2 + I * K * ( 2/(step_size/2) - n[ik].imag() );
          /*
           * absorption borders are formed here 
           */
          if ( i <= i_left )
            ci[i].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(i_left - i) /
                                     static_cast<double>(i_left),
                                     band_pow );

          if ( i >= i_right )
            ci[i].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(i - i_right) /
                                     static_cast<double>(info.number.first - i_right),
                                     band_pow );
        }

        for ( size_t i = 0; i < info.number.first; ++i ) {
          idx(i,j-1)  = ui[i];
        }

        elim_i(ci, pi, ui);

        j = jj + 1;
        for ( size_t i = 1; i < info.number.first - 1; ++i ) {
          const Pixel & uij   = idx(i,j    );
          const Pixel & uij1  = idx(i,j + 1);
          const Pixel & uij_1 = idx(i,j - 1);

          pi[i-1] = 2. * I * K / (step_size/2.0) * uij
                  - (-2. * uij + uij1 + uij_1) / dy2;
        }

        for ( size_t i = 0; i < info.number.first; ++i ) {
          const size_t ik = i * info.number.second + j;

          ci[i] = -2/dy2 + I * K * ( 2/(step_size/2) - n[ik].imag() );
          /*
           * absorption borders are formed here 
           */
          if ( i <= i_left )
            ci[i].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(i_left - i) /
                                     static_cast<double>(i_left),
                                     band_pow );

          if ( i >= i_right )
            ci[i].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(i - i_right) /
                                     static_cast<double>(info.number.first - i_right),
                                     band_pow );
        }

        for ( int i = 0; i < info.number.first; ++i ) {
          idx(i,j-1) = ui[i];
        }

        elim_i(ci, pi, ui);
      }

      /* something to the last column. */
      for ( int i = 0; i < info.number.first; ++i ) {
        idx(i,info.number.second-1) = ui[i];
      }
      

      /*
       * end j 
       */

      /***************************************************/


      if ( dump_file.is_open() &&
           ( (istep/dump_period) ==
             static_cast<float>(istep)/static_cast<float>(dump_period)
           )
         ) {
        /* writing out a cross-section of the optical field.  This is
         * essentially a copy of cros_out.cpp */

        /*
         * writing the intensity of x-direction
         */
        const double dx = info.side_length.first / (info.number.first - 1.);
        const double dy = info.side_length.second / (info.number.second - 1.);
        const size_t i2 = info.number.first / 2;
        const size_t j2 = info.number.first / 2;
        const double z = istep * step_size;

        for ( size_t i = 0; i < info.number.first; ++i ) {
          double x = dx * ( i - (double)i2 );
          const Pixel & px  = idx(i,j2);
          dump_file << z << '\t' << x << '\t'
                    << norm(px) << '\t' << arg(px) << '\n';
        }

        dump_file << '\n'; // blank line between x and y directions

        /*
         * writing the intensity of y-direction
         */
        for ( int j = 0; j < info.number.second; ++j ) {
          double y = dy * ( j - (double)j2 );
          const Pixel & px  = idx(i2,j);
          dump_file << z << '\t' << y << '\t'
                    << norm(px) << '\t' << arg(px) << '\n';
        }

        dump_file << std::endl;
      }
    }

    for ( int i = info.size() - 1; i >= 0; --i ) {
      register double fi = 0.25 * K * (step_size/2.0) * ( n[i].real() - 1. );
      val[i] *= exp(I * fi);
    }

    if ( dump_file.is_open() ) {
      dump_file.close();
    }
  }
  /* *** STUFF FOR STEPS ENDS HERE *** */


  void Field::cleanup() {
    if (val) {
      delete[] val;
      val = NULL;

      fftw_destroy_plan(fftw_configs[0]);
      fftw_destroy_plan(fftw_configs[1]);
    }
  }

  void Field::init(Field::Pixel init_fill) {
    cleanup();
    this->base_init();
    std::fill( val, val + info.size(), init_fill );
  }

}/* namespace lightpipes */
