

#include <lightpipes/Field.h>
extern "C" {
  #include <lightpipes/fftn.h>
}

#include <algorithm>
#include <fstream>

#include <stdint.h>
#include <cmath>


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


  Field::Field ( const std::pair<size_t,size_t> & number,
                 const std::pair<double,double> & side_length,
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


  Field & Field::operator= (const Field & that) {
    cleanup();
    this->info = that.info;
    size_t sz = info.size();

    val = new std::complex<double>[sz];
    if (val == NULL)
      throw std::runtime_error("field allocation error");

    std::copy( that.val, that.val+sz, val );
    return *this;
  }


  Field & Field::operator+= ( const Field & that ) {
    if (!compatible(that))
      throw std::runtime_error("cannot add fields that are not compatible");

    for ( int i = info.size()-1; i >= 0; --i ) {
      val[i] += that.val[i];
    }

    return *this;
  }



  Field & Field::circular_aperture(const double & r, const double & x0, const double & y0) {
    double r2 = SQR(r);
    double dx = info.side_length.first / info.number.first;
    double dy = info.side_length.second / info.number.second;
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    /*
     * Cutting the aperture 
     */

    for ( int i = 0; i < info.number.first; ++i ) {
      double x2  = SQR((i - i2) * dx - x0);
      for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    for ( int i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);
      for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    for ( int i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);
      for ( int j = 0; j < info.number.second; ++j ) {
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
      for ( int i = sz - 1; i >= 0; --i )
        val[i] /= ampl_scale;
    } else {
      Pixel * f1 = new Pixel[sz];
      if ( f1 == NULL ) {
        throw std::runtime_error("Allocation error, f1 in lens_forvard");
      }

      std::fill(f1, f1+sz, Pixel(0.0));

      for ( int i = 0; i < info.number.first; ++i ) {
        size_t i_i = info.number.first - i + 1;
        for ( int j = 0; j < info.number.second; ++j ) {
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
      int i2 = info.number.first / 2;
      int j2 = info.number.second / 2;
      int sign = (z >= 0.) ? 1 : -1;

      for ( int i = 0; i < info.number.first; ++i ) {
        for ( int j = 0; j < info.number.second; ++j ) {
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
    int dims[2] = { info.number.first, info.number.second };

    size_t sz = info.size();
    double * Re = new double[sz];
    double * Im = new double[sz];
    for (int i = 0; i < sz; ++i) {
      Re[i] = val[i].real();
      Im[i] = val[i].imag();
    }

    /** FIXME:  Not sure if this scaling was correct
     * the scaling, before switching to NxM was (double)info.number */
    fftn ( 2, dims, Re, Im, direction, -2 );

    for (int i = 0; i < sz; ++i)
      val[i] = Pixel (Re[i], Im[i]);

    delete[] Re;
    delete[] Im;
    return *this;
  }

  Field & Field::circular_screen( const double & r, const double & x0, const double & y0 ) {
    double r2 = SQR(r);

    double dx = info.side_length.first / ( info.number.first );
    double dy = info.side_length.second / ( info.number.second );
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    /*
     * Cutting the aperture 
     */

    for ( int i = 0; i < info.number.first; ++i ) {
        double x = (i - i2) * dx - x0;
        for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    double cc = cos ( angle );
    double ss = sin ( angle );

    /*
     * Cutting the aperture 
     */

    for ( int i = 0; i < info.number.first; ++i ) {
      double xx = (i - i2) * dx - x0;
      for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    double cc = cos ( angle );
    double ss = sin ( angle );

    /*
     * Cutting the aperture 
     */

    for ( int i = 0; i < info.number.first; ++i ) {
      double xx = (i - i2) * dx - x0;
      for ( int j = 0; j < info.number.second; ++j ) {
        double yy = (j - i2) * dy - y0;

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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    double w2_inv = 1./ (SQR(w) * 2.);
    double Asqrt = std:sqrt ( abs(A) );

    for ( int i = 0; i < info.number.first; ++i ) {
        double x2 = SQR( (i - i2) * dx - x0 );
        for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    double w2_inv = 1./SQR(w);

    for ( int i = 0; i < info.number.first; ++i ) {
        double x2 = SQR( (i - i2) * dx - x0 );
        for ( int j = 0; j < info.number.second; ++j ) {
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
    for ( int i = 0; i < info.number.first; ++i ) {
      for ( int j = 0; j < info.number.second; ++j ) {
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
                          const std::pair<double,double> & new_side_length,
                          const std::pair<size_t,size_t> & new_number ) {
    throw std::runtime_error("not fixed for NxM yet");
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
  }

  Field & Field::fresnel ( const double &z ) {
    throw std::runtime_error("not fixed for NxM yet");
    int ii, ij;

    /*
     * fprintf(stderr,"%e %d %e\n", z,n_old, dx); 
     */
    double dx = info.side_length / ( info.number - 1. );
    /*
     * Allocating a LOT OF MEMORY 
     */

    int fn2 = info.number * 2;

    int len = SQR(fn2);

    double * F_R = new double[len];
    if ( F_R == NULL ) {
        throw std:: runtime_error ( "Allocation error in Fresnel, use smaller grid" );
    }

    double * F_I = new double[len];
    if ( F_I == NULL ) {
        throw std:: runtime_error ( "Allocation error in Fresnel, use smaller grid" );
    }

    double * K_R = new double[len];
    if ( K_R == NULL ) {
        throw std:: runtime_error ( "Allocation error in Fresnel, use smaller grid" );
    }

    double * K_I = new double[len];
    if ( K_I == NULL ) {
        throw std:: runtime_error ( "Allocation error in Fresnel, use smaller grid" );
    }

    std::fill( F_R, F_R+len, 0.0 );
    std::fill( F_I, F_I+len, 0.0 );
    std::fill( K_R, K_R+len, 0.0 );
    std::fill( K_I, K_I+len, 0.0 );

    /*****************************************************/

    double sh = 0.5;
    int fn22 = info.number + 1;
    int no2 = info.number / 2;
    double RR = sqrt ( 1. / ( 2 * info.lambda * z ) ) * dx * 2;

    ii = ij = 1;
    long ik = 0;

    for ( int i = fn22 - no2; i <= fn22 + no2 - 1; i++ ) {
        int io = i - fn22;

        double R1 = RR * ( io - .5 + sh );
        double R3 = RR * ( io + .5 + sh );

        for ( int j = fn22 - no2; j <= fn22 + no2 - 1; j++ ) {
            int iiij = ii * ij;
            int jo = j - fn22;
            long ik1 = ( i - 1 ) * fn2 + j - 1;

            /*
             * Fresnel staff 
             */
            double R2 = RR * ( jo - .5 + sh );
            double R4 = RR * ( jo + .5 + sh );

            double fc1(0.), fs1(0.),
                   fc2(0.), fs2(0.),
                   fc3(0.), fs3(0.),
                   fc4(0.), fs4(0.);

            fresnl ( R1, &fs1, &fc1 );
            fresnl ( R2, &fs2, &fc2 );
            fresnl ( R3, &fs3, &fc3 );
            fresnl ( R4, &fs4, &fc4 );

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


            K_R[ik1] = 0.5 * ( c4s3 + s4c3 - c4s1 - s4c1 - c2s3 - s2c3 + c2s1 + s2c1 ) * iiij;
            K_I[ik1] = 0.5 * ( -c4c3 + s4s3 + c4c1 - s4s1 + c2c3 - s2s3 - c2c1 + s2s1 ) * iiij;
            /*
             * Field staff 
             */

            /*
             * F_R[ik1]=(field.real[ik]-field.real[ik-1]+field.real[ik-info.number]-field.real[ik-info.number-1])*iiij*0.25;
             * F_I[ik1]=(field.imaginary[ik]-field.imaginary[ik-1]+field.imaginary[ik-info.number] -
             * field.imaginary[ik-info.number-1])*0.25*iiij; 
             */

            F_R[ik1] = ( val[ik].real() ) * iiij;
            F_I[ik1] = ( val[ik].imag() ) * iiij;

            ik++;
            ij = -ij;
        }
        ii = -ii;
    }

    int dims[2] = {fn2, fn2};

    fftn ( 2, dims, K_R, K_I, 1, ( double ) fn2 );
    fftn ( 2, dims, F_R, F_I, 1, ( double ) fn2 );

    ik = 0;
    ii = ij = 1;
    for ( int i = 1; i <= fn2; i++ ) {

        for ( int j = 1; j <= fn2; j++ ) {
            int iiij = ii * ij;

            double cc = K_R[ik] * F_R[ik] - K_I[ik] * F_I[ik];
            F_I[ik] = ( K_R[ik] * F_I[ik] + F_R[ik] * K_I[ik] ) * iiij;
            F_R[ik] = cc * iiij;
            ik++;
            ij = -ij;
        }
        ii = -ii;
    }

    delete[] K_R;
    delete[] K_I;

    fftn ( 2, dims, F_R, F_I, -1, 1. );

    ik = 0;
    ii = ij = 1;
    for ( int i = fn22 - no2; i <= fn22 + no2 - 1; i++ ) {

        for ( int j = fn22 - no2; j <= fn22 + no2 - 1; j++ ) {
            long ik1 = ( i - 1 ) * fn2 + j - 1;
            long ik2 = ( i - 2 ) * fn2 + j - 1;
            long ik3 = ( i - 2 ) * fn2 + j - 2;
            long ik4 = ( i - 1 ) * fn2 + j - 2;
            int iiij = ii * ij;

            double Fr = 0.25 * ( F_R[ik1] - F_R[ik2] + F_R[ik3] - F_R[ik4] ) * iiij;
            double Fi = 0.25 * ( F_I[ik1] - F_I[ik2] + F_I[ik3] - F_I[ik4] ) * iiij;

            val[ik] = std::complex<double>(Fr,Fi);

            /*
             * field.real[ik]=F_R[ik1]*iiij;
             * field.imaginary[ik]=F_I[ik1]*iiij; 
             */
            ik++;
            ij = -ij;
        }
        ii = -ii;
    }

    delete[] F_R;
    delete[] F_I;

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
    int i2 = info.number.second / 2;
    int j2 = info.number.second / 2;

    for ( int i = 0; i < info.number.first; ++i ) {
      double x = (i - i2) * dx;
      for ( int j = 0; j < info.number.second; ++j ) {
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
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    for ( int i = 0; i < info.number.first; ++i ) {
      double x = (i - i2) * dx;
      for ( int j = 0; j < info.number.second; ++j ) {
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
    for ( int i = 0; i < info.number.first; ++i ) {
      for ( int j = 0; j < info.number.second; ++j ) {
        Pixel & px = idx(i,j);
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
    double dx2 = SQR(dx);
    double dy2 = SQR(dy);
    int i2 = info.number.first / 2;
    int j2 = info.number.second / 2;

    { std::pair<double,double> strehl_energy = get_strehl_and_energy();
      output << "Strehl: ratio= " << strehl_energy.first
             << " energy= " << (strehl_energy.second * dx2)
             << std::endl;
    }

    /*
     * Calculating the center of gravity: 
     */
    double sum = 0., sum1r = 0., sum1i = 0.;
    for ( int i = 0; i < info.number.first; ++i ) {
      double y = (i - i2) * dx;
      for ( int j = 0; j < info.number.second; ++j ) {
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
    for ( int i = 0; i < info.number.first; ++i ) {
        double y = (i - i2) * dx;
        double y_y_c = y - y_c;
        for ( int j = 0; j < info.number.second; ++j ) {
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
                                    Field::SzPair   output_size /* = 0 */,
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

    int istep = 1, jstep = 1;

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
      for ( int i = info.size()-1; i >= 0; --i ) {
        double sum = norm( val[i] );
        if ( sum > max_px ) max_px = sum;
      }
    }

    /* now write pixel (set) values out to file. */
    int i_i = 1;
    for ( int i = 0; i <= info.number.first - istep; i += istep ) {
      for ( int j = 0; j <= info.number.second - jstep; j += jstep ) {
        /* determine the average intensity/phase in the next pixel set. */

        double aval = 0.0;
        if (output_norm) {
          double tot = 0;
          for ( int ii = i; ii < i + istep; ++ii )
            for ( int jj = j; jj < j + jstep; ++jj )
              tot += norm( idx(ii,jj) );

          aval = tot / (istep * jstep);
        } else {
          Pixel tot = 0;
          for ( int ii = i; ii < i + istep; ++ii )
            for ( int jj = j; jj < j + jstep; ++jj )
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
    for ( int i = info.size()-1; i >= 0; --i )
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
    for ( int i = info.size()-1; i >= 0; --i )
      val[i] *= sum;

    return *this;
  }

  Field & Field::l_amplify(const double & gain,
                           const double & length,
                           const double & i_sat) {
    using std::norm; using std::exp;
    for ( int i = info.size()-1; i >= 0; --i ) {
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
    double dx    = info.side_length.first  / info.number.first;
    double dy    = info.side_length.second / info.number.second;
    int    i2    = info.number.first  / 2;
    int    j2    = info.number.second / 2;

    for ( int i = 0; i < info.number.first; ++i ) {
      double x2 = SQR((i - i2) * dx - x0);

      for ( int j = 0; j < info.number.second; ++j ) {
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



  Field & Field::interpolate(const std::pair<double,double> & new_side_length /* = 0.0 */,
                             const SzPair & new_number /* = 0.0 */,
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
  }

  /* *** STUFF FOR INTERPOLATE ENDS HERE *** */

  namespace {
    void elim ( const Field::Info & info,
                      std::complex<double> * alpha,
                      std::complex<double> * beta,
                const std::complex<double> * a,
                const std::complex<double> * b,
                const std::complex<double> * c,
                const std::complex<double> * p,
                      std::complex<double> * u ) {
      throw std::runtime_error("not fixed for NxM yet");
      /*
       * initial condition, everything is going to be zero at the edge 
       */
      alpha[info.number] = beta[info.number]
        = alpha[2] = beta[2] = std::complex<double>(0.0, 0.0);

      /*
       * forward elimination 
       */
      std::complex<double> cc;
      for ( int i = 2; i <= info.number - 2; ++i ) {
        cc = c[i] -  a[i] * alpha[i];
        alpha[i + 1] = b[i] / cc;
        beta[i + 1] = ( p[i] + a[i] * beta[i] ) / cc;
      }
      cc = c[info.number] - a[info.number] * alpha[info.number];
      beta[info.number + 1]
        = ( p[info.number] + ( a[info.number] * beta[info.number] ) ) / cc;
      /*
       * edge amplitude =0 
       */
      u[info.number] = beta[info.number + 1];
      /*
       * backward elimination 
       */
      for ( int i = info.number - 1; i >= 1; --i ) {
        u[i] = ( alpha[i + 1] * u[i + 1] ) + beta[i + 1];
      }
    }
  }

  /* *** STUFF FOR STEPS BEGINS HERE *** */
  Field & Field::steps( const double & step_size,
                        const int & N,
                        const std::string & n_filename,
                        const std::string & k_filename,
                        const std::string & dump_filename,
                        const int & dump_period)
                        throw (std::runtime_error) {
    throw std::runtime_error("not fixed for NxM yet");
    typedef std::complex<double> Complex;

    Complex * n = new Complex[ SQR(info.number) ];
    if ( n == NULL )
      throw std::runtime_error("Field::step: Allocation error for 'n'");
    std::fill( n, n+SQR(info.number), Complex(1.0, 0.0) );

    if ( n_filename != "" ) {
      std::ifstream n_file( n_filename.c_str() );
      if ( !n_file )
        throw std::runtime_error(
          "Field::steps: error opening file refractive index file (%s)"
        );

      long ik = 0;
      for ( int i = 1; i <= info.number; ++i ) {
        for ( int j = 1; j <= info.number; ++j ) {
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


      long ik = 0;
      for ( int i = 1; i <= info.number; ++i ) {
        for ( int j = 1; j <= info.number; ++j ) {
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

  Field & Field::steps( const double & step_size,
                        const int & N,
                        const std::complex<double> * n,
                        const std::string & dump_filename,
                        const int & dump_period)
                        throw (std::runtime_error) {
    throw std::runtime_error("not fixed for NxM yet");
    typedef std::complex<double> Complex;
    const double K = 2. * M_PI / info.lambda;

    bool clear_n = false;
    if ( n == NULL ) {
      clear_n = true;
      Complex * tmp = new Complex[ SQR(info.number) ];
      if ( tmp == NULL )
        throw std::runtime_error("Field::step: Allocation error for 'n'");
      std::fill( tmp, tmp+SQR(info.number), Complex(1.0, 0.0) );
      n = tmp;
    }

    Complex * alpha= new Complex[ info.number+3 ];
    Complex * beta = new Complex[ info.number+3 ];
    Complex * a    = new Complex[ info.number+3 ];
    Complex * b    = new Complex[ info.number+3 ];
    Complex * c    = new Complex[ info.number+3 ];
    Complex * p    = new Complex[ info.number+3 ];
    Complex * u    = new Complex[ info.number+3 ];
    Complex * u1   = new Complex[ info.number+3 ];
    Complex * u2   = new Complex[ info.number+3 ];


    std::ofstream dump_file;
    double * int1 = NULL;
    double * phase1 = NULL;
    if ( dump_filename != "" ) {
      dump_file.open( dump_filename.c_str() );

      int1 = new double[info.number + 2];
      if ( int1 == NULL )
        throw std::runtime_error("Field::steps: Allocation error for int1");

      std::fill( int1, int1+(info.number + 2), 0.0 );

      phase1 = new double[ info.number + 2 ];
      if ( phase1 == NULL )
        throw std::runtime_error("Field::steps: Allocation error for phase1");

      std::fill( phase1, phase1+(info.number + 2), 0.0 );
    }

    /*
     * the arrays for refraction and absorption are finally formed 
     */

    const double step_size_2 = step_size / 2.;
    const double Pi4lz = 4. * M_PI / info.lambda / step_size_2;
    const double delta2 = SQR( info.side_length / (info.number - 1.0 ) );
    /*
     * absorption at the borders is described here 
     */
    const double AA = -10. / step_size_2 / N; /* total absorption */
    const double band_pow = 2.; /* profile of the absorption border,
                                 * 2=quadratic */
    /*
     * width of the absorption border 
     */
    const int i_left = info.number / 2 + 1 - 0.4 * info.number;
    const int i_right = info.number / 2 + 1 + 0.4 * info.number;
    /*
     * end absorption 
     */


    for ( int i = 1; i <= info.number; ++i ) {
      u2[i] = Complex(0.0, 0.0);
      a[i] = b[i] = Complex(-1.0 / delta2, 0.0);
    }

    double dist = 0.;
    /*
     * Main loop, steps here 
     */
    for ( int istep = 1; istep <= N; istep++ ) {
      dist += 2. * step_size_2;
      /*
       * Elimination in the direction i, halfstep 
       */

      { int ik = 0;
        for ( int i = 1; i <= info.number; ++i ) {
          for ( int j = 1; j <= info.number; ++j ) {
            register double fi = 0.25 * K * step_size_2 * ( n[ik].real() - 1. );
            val[ik++] *= exp(I * fi);
          }
        }
      }


      for ( int jj = 2; jj <= info.number - 2; jj += 2 ) {
        int j = jj;

        for ( int i = 2; i <= info.number - 1; ++i ) {
          const int ikij   = ( j - 1 ) * info.number + i - 1;
          const int ikij1  = ( j     ) * info.number + i - 1;
          const int ikij_1 = ( j - 2 ) * info.number + i - 1;

          const Complex & uij   = val[ikij];
          const Complex & uij1  = val[ikij1];
          const Complex & uij_1 = val[ikij_1];

          p[i] = -2. * uij;
          p[i] = p[i] + uij1;
          p[i] = p[i] + uij_1;
          p[i] = ( -1. / delta2 ) * p[i];
          p[i] =  Complex ( 0., Pi4lz ) * uij +  p[i];
        }

        for ( int i = 1; i <= info.number; i++ ) {
          const int ik = ( j - 1 ) * info.number + i - 1;

          const double medium_i = -2. * M_PI * n[ik].imag() / info.lambda;
          c[i].real() = -2. / delta2;
          c[i].imag() = Pi4lz + medium_i;
          /*
           * absorption borders are formed here 
           */
          if ( i <= i_left )
            c[i].imag() -= ( AA * K )
                         * std::pow ( static_cast<double>(i_left - i) /
                                      static_cast<double >(i_left),
                                      band_pow );

          if ( i >= i_right )
            c[i].imag() -= ( AA * K )
                         * std::pow( static_cast<double>(i - i_right) /
                                     static_cast<double>(info.number - i_right),
                                     band_pow );
          /*
           * end absorption 
           */

        }


        elim(info, alpha, beta, a, b, c, p, u);


        for ( int i = 1; i <= info.number; ++i ) {
          const int ikij_1 = ( j - 2 ) * info.number + i - 1;
          val[ikij_1] = u2[i];
          u1[i] = u[i];
        }

        j = jj + 1;

        for ( int i = 2; i <= info.number - 1; ++i ) {
          const int ikij   = ( j - 1 ) * info.number + i - 1;
          const int ikij1  = ( j     ) * info.number + i - 1;
          const int ikij_1 = ( j - 2 ) * info.number + i - 1;

          const Complex & uij   = val[ikij];
          const Complex & uij1  = val[ikij1];
          const Complex & uij_1 = val[ikij_1];

          p[i] = -2. * uij;
          p[i] = p[i] + uij1;
          p[i] = p[i] + uij_1;
          p[i] = ( -1. / delta2) * p[i];
          p[i] = Complex ( 0., Pi4lz ) * uij + p[i];
        }

        for ( int i = 1; i <= info.number; ++i ) {
          const int ik = ( j - 1 ) * info.number + i - 1;

          const double medium_i = -2. * M_PI * n[ik].imag() / info.lambda;
          c[i].real() = -2. / delta2;
          c[i].imag() = Pi4lz + medium_i;

          /*
           * absorption borders are formed here 
           */
          if ( i <= i_left )
            c[i].imag() -= ( AA * K )
                         * std::pow ( static_cast<double>(i_left - i) /
                                      static_cast<double>(i_left),
                                      band_pow );

          if ( i >= i_right )
            c[i].imag() -= ( AA * 2. * K )
                         * std::pow ( static_cast<double>(i - i_right) /
                                      static_cast<double>(info.number - i_right),
                                      band_pow );
          /*
           * end absorption 
           */

        }


        elim(info, alpha, beta, a, b, c, p, u);


        for ( int i = 1; i <= info.number; ++i ) {
          const int ikij = ( j - 2 ) * info.number + i - 1;
          val[ikij] = u1[i];
          u2[i] = u[i];
        }
      }

      for ( int i = 1; i <= info.number; i++ ) {
        const int ikij = ( info.number - 1 ) * info.number + i - 1;
        val[ikij] = u2[i];
      }

      { int ik = 0;
        for ( int i = 1; i <= info.number; ++i ) {
          for ( int j = 1; j <= info.number; ++j ) {
            register double fi = 0.5 * K * step_size_2 * ( n[ik].real() - 1. );
            val[ik++] *= exp(I * fi);
          }
        }
      }

      /*
       * Elimination in the j direction is here, halfstep 
       */

      for ( int i = 1; i <= info.number; ++i ) {
        u2[i] = Complex(0.0, 0.0);
      }

      { int i = 0;
        for ( int ii = 2; ii <= info.number - 2; ii += 2 ) {
          i = ii;
          for ( int j = 2; j <= info.number - 1; ++j ) {
            const int ikij   = ( j - 1 ) * info.number + i - 1;
            const int iki1j  = ( j - 1 ) * info.number + i;
            const int iki_1j = ( j - 1 ) * info.number + i - 2;

            const Complex & uij   = val[ikij];
            const Complex & ui1j  = val[iki1j];
            const Complex & ui_1j = val[iki_1j];

            p[j] = -2. * uij;
            p[j] = p[j] + ui1j;
            p[j] = p[j] + ui_1j;
            p[j] = ( -1. / delta2 ) * p[j];
            p[j] = Complex ( 0., Pi4lz ) * uij + p[j];
          }

          for ( int j = 1; j <= info.number; ++j ) {
            const int ik = ( j - 1 ) * info.number + i - 1;

            const double medium_i = -2. * M_PI * n[ik].imag() / info.lambda;
            c[j].real() = -2. / delta2;
            c[j].imag() = Pi4lz + medium_i;

            /*
             * absorption borders are formed here 
             */
            if ( j <= i_left )
              c[j].imag() -= ( AA * K )
                           * std::pow ( static_cast<double>(i_left - j) /
                                        static_cast<double>(i_left),
                                        band_pow );

            if ( j >= i_right )
              c[j].imag() -= ( AA * K )
                           * std::pow ( static_cast<double>(j - i_right) /
                                        static_cast<double>(info.number - i_right),
                                        band_pow );

            /*
             * end absorption 
             */

          }


          elim(info, alpha, beta, a, b, c, p, u);


          for ( int j = 1; j <= info.number; ++j ) {
            const int iki_1j = ( j - 1 ) * info.number + i - 2;
            val[iki_1j] = u2[j];
            u1[j] = u[j];
          }

          i = ii + 1;
          for ( int j = 2; j <= info.number - 1; ++j ) {
            const int ikij   = ( j - 1 ) * info.number + i - 1;
            const int iki1j  = ( j - 1 ) * info.number + i;
            const int iki_1j = ( j - 1 ) * info.number + i - 2;

            const Complex & uij = val[ikij];
            const Complex & ui1j = val[iki1j];
            const Complex & ui_1j = val[iki_1j];

            p[j] = -2. * uij;
            p[j] = p[j] + ui1j;
            p[j] = p[j] + ui_1j;
            p[j] = ( -1. / delta2 ) * p[j];
            p[j] = Complex ( 0., Pi4lz ) * uij + p[j];
          }

          for ( int j = 1; j <= info.number; ++j ) {
            const int ik = ( j - 1 ) * info.number + i - 1;

            const double medium_i = -2. * M_PI * n[ik].imag() / info.lambda;
            c[j].real() = -2. / delta2;
            c[j].imag() = Pi4lz + medium_i;

            /*
             * absorption borders are formed here 
             */
            if ( j <= i_left )
              c[j].imag() -= ( AA * K )
                           * std::pow ( static_cast<double>(i_left - j) /
                                        static_cast<double>(i_left),
                                        band_pow );

            if ( j >= i_right )
              c[j].imag() -= ( AA * K )
                           * std::pow ( static_cast<double>(j - i_right) /
                                        static_cast<double>(info.number - i_right),
                                        band_pow );

            /*
             * end absorption 
             */

          }


          elim(info, alpha, beta, a, b, c, p, u);


          for ( int j = 1; j <= info.number; ++j ) {
            const int ikij = ( j - 1 ) * info.number + i - 2;
            val[ikij] = u1[j];
            u2[j] = u[j];
          }
        }

        for ( int j = 1; j <= info.number; ++j ) {
          const int ikij = ( j - 1 ) * info.number + i - 1;
          val[ikij] = u2[j];
        }
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

        /*
         * writing the intensity into arrays 
         */

        { const int i = info.number / 2 + 1;
          for ( int j = 1; j <= info.number; j += 1 ) {
            const int ik1 = ( i - 1 ) * info.number + j - 1;
            const int jj = j - 1;
            int1[jj] = norm( val[ik1] );
            phase1[jj] = arg(val[ik1]);
          }
        }

        { const int j = info.number / 2 + 1;
          const double dx = info.side_length / ( info.number - 1. );
          for ( int i = 1; i <= info.number; i += 1 ) {
            const int ik1 = ( i - 1 ) * info.number + j - 1;
            const int jj = i - 1;
            double cc = dx * ( i - info.number / 2 - 1 );

            double int2   = norm( val[ik1] );
            double phase2 =  arg( val[ik1] );
            dump_file << cc << '\t'
                      << int1[jj] << '\t'
                      << int2 << '\t'
                      << phase1[jj] << '\t'
                      << phase2 << '\t'
                      << dist << '\n';
          }
        }

        dump_file << std::endl;
      }
    }


    { int ik = 0;
      for ( int i = 1; i <= info.number; ++i ) {
        for ( int j = 1; j <= info.number; ++j ) {
          register double fi = 0.25 * K * step_size_2 * ( n[ik].real() - 1. );
          val[ik++] *= exp(I * fi);
        }
      }
    }

    if (clear_n)
      delete[] n;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] u;
    delete[] u1;
    delete[] u2;
    delete[] alpha;
    delete[] beta;
    delete[] p;

    if ( dump_file.is_open() ) {
      delete[] int1;
      delete[] phase1;
    }

    return *this;
  }
  /* *** STUFF FOR STEPS ENDS HERE *** */


  void Field::cleanup() {
    if (val) {
      delete[] val;
      val = NULL;
    }
  }

  void Field::init(Field::Pixel init_fill) {
    cleanup();

    val = new Pixel[info.size()](init_fill);
    if (val == NULL)
      throw std::runtime_error("field allocation error");
  }

}/* namespace lightpipes */
