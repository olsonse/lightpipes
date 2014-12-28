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

#ifndef lightpipies_Field_h
#define lightpipies_Field_h

#include <fftw3.h>

#include <iostream>
#include <stdexcept>
#include <complex>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef _DJGPP_
#  include <fcntl.h>
#endif

namespace lightpipes {

  template < typename T> T SQR(const T & t) { return t*t; }


  template < typename T1, typename T2=T1 >
  struct Pair : std::pair<T1,T2> {
    typedef std::pair<T1,T2> super;
    Pair( const std::pair<T1,T2> & that ) {
      (*this) = static_cast<const Pair&>(that);
    }
    Pair( const Pair & that ) {
      (*this) = that;
    }
    Pair( const T1 & i, const T2 & j ) : super(i,j) { }
    Pair( const T1 & i = T1(0) ) : super(i,i) { }
  };

  template < typename T1, typename T2 >
  std::ostream & operator<<( std::ostream & out, const Pair<T1,T2> & p ) {
    return out << '(' << p.first << ',' << p.second << ')';
  }


#ifdef _DJGPP_
  extern void setmode ( int, int );
#endif


  /** Field structure represents the light field at a given place in
   * time/position.
   *
   * The structure FIELD contains the characteristics of the light beam:
   * number of points along the side of a square grid, wavelength and side
   * length of the square grid, then two huge arrays of Re and Im data 
   */
  class Field {
    /* TYPEDEFS */
  public:
    /** Type of pixel in Field. */
    typedef std::complex< double > Pixel;

    /** Defines character of Field such as wavelength of light, spatial size,
     * spatial precision, etc. */
    struct Info {
      /* STORAGE MEMBERS */
      /** The number of side-elements stored in Field::val
       * (sizeof(Field::val) = number.first * number.second.<br>
       * number.first denotes the number of rows<br>
       * number.second denotes the number of columns
       */
      Pair<size_t> number;

      /** The physical size of the side of the Field view. */
      Pair<double> side_length;

      /** Wavelength of the Field. */
      double lambda;

      /** FFT status. */
      int fft_level;

      /** Spherical coordinates conversion factor, see Field::lens_fresnel
       * and Field::lens_forvard. */
      double sph_coords_factor;



      /* FUNCTION MEMBERS */
      /** Constructor--zeros everything. */
      Info();

      /** Test to see if two Fields are compatible based on their field type
       * (wavelength, size, etc.).
       */
      bool compatible(const Info & that) const;

      /** The number of pixels in the field. */
      inline size_t size() const {
        return number.first * number.second;
      }

      /** The number of pixels in the field. */
      inline size_t mem_size() const {
        return this->size() * sizeof(Pixel);
      }

      /** The number of pixels in the field. */
      inline double phy_size() const {
        return side_length.first * side_length.second;
      }

      /** read the Field::Info data from file. */
      std::istream & read(std::istream & in);

      /** write the Field::Info data to file. */
      std::ostream & write(std::ostream & out);
    };



    /* STORAGE MEMBERS */
  public:
    /** Field::Info describes wavelength, view size, view area. */
    Info info;

    /** Array of Field values. */
    Pixel * val;

    /** FFTW3 Plan */
    fftw_plan fftw_configs[2];



    /* FUNCTION MEMBERS */
  public:
    /** Constructor. */
    Field ( const Pair<size_t> & number,
            const Pair<double> & side_length,
            double lambda,
            Pixel init_fill = Pixel(1.,0.),
            int fft_level = 0,
            double sph_coords_factor = 0.0 );

    /** Constructor with Field::Info initialization. */
    Field (const Info & that_info,
           Pixel init_fill = Pixel(1.,0.) );

    /** Copy constructor. */
    Field (const Field & that);

    /** Destructor. */
    ~Field ();

    /** Copy generator. */
    Field copy() const {
      return Field(*this);
    }

    /** Write the Field to file (internal format). 
    * The internal format is 1. Field::Info, 2. complex<double> array of
    * length N*N.
    */
    std::ostream & write(std::ostream & out = std::cout);

    /** Read a new Field from file (internal format). 
     * @see Field::write.
     */
    static Field * read(std::istream & in = std::cin) throw (std::runtime_error) ;

    /** Field index operator (non-const). */
    Pixel & operator[](const size_t i) { return val[i]; }

    /** Field index operator (const). */
    const Pixel & operator[](const size_t i) const { return val[i]; }

  private:
    /** Field index operator (non-const) using column AND row indices. */
    Pixel & idx(const size_t row, const size_t col) {
      return val[row*info.number.second + col];
    }

    /** Field index operator (const) using column AND row indices. */
    const Pixel & idx(const size_t row, const size_t col) const {
      return val[row*info.number.second + col];
    }


  public:
    /** Field index operator (non-const) using column AND row indices. */
    Pixel & operator()(const size_t row, const size_t col) {
      return idx(row,col);
    }

    /** Field index operator (const) using column AND row indices. */
    const Pixel & operator()(const size_t row, const size_t col) const {
      return idx(row,col);
    }

    /** Copy operator. */
    Field & operator= (const Field & that);

    /** fill the field with a constant value operator. */
    template < typename T >
    Field & fill(const T & t) {
      std::fill( val, val + info.size(), t );
      return *this;
    }

    /** Multiply alternating field elements by (-1)^(row) */
    Field & negate_alternate_elems() {
      int ii = 1, ij = 1;
      for ( size_t i = 0; i < info.number.first; ++i ) {
        for ( size_t j = 0; j < info.number.second; ++j ) {
          idx(i,j) *= ii * ij;
          ij = -ij;
        }
        ii = -ii;
      }
      return *this;
    }

    /** Scaling operator. */
    template <class T>
    Field & operator*= ( const T & input ) {
      for (long i = info.size()-1; i >= 0; --i)
        val[i] *= input;
      return *this;
    }

    /** Field addition operator (Field1 += Field2). */
    Field & operator+= ( const Field & that );


    /** Apply a lens.
     * f>0 for positive lens !!!!!  (Unlike the original lightpipes code)
     */
    Field & lens ( const double & f,
                   const double & x0 = 0.0,
                   const double & y0 = 0.0 );

    /** Apply a toroidal lens.
     * f>0 for positive lens !!!!!  (Unlike the original lightpipes code)
     * filters the beam through the toroidal quadratic phase corrector.
     *
     * @param R toroidal radius.
     * @param f focal length.
     * @param x0 x-shift in position of center.
     * @param y0 y-shift in position of center.
     */
    Field & t_lens ( const double & R,
                     const double & f,
                     const double & x0 = 0.0,
                     const double & y0 = 0.0 );

    /** Apply an axicon lens.
     * @param phi including angle of axicon (note that the sign of this value
     * is disregarded).
     * @param n1 complex index of refraction of axicon material.
     * @param x0 x-shift in position of center.
     * @param y0 y-shift in position of center.
     */
    Field & axicon ( const double & phi,
                     const Pixel & n1 = Pixel(1.5,0.0),
                     const double & x0 = 0.0,
                     const double & y0 = 0.0 );

    /** Propagates Field using convolution.
     * Another possibility of a fast computer implementation of the operator
     * \f$L^+\f$ is free from many of the drawbacks of the
     * described spectral algorithm. The operator \f$L^+\f$ may
     * be numerically implemented with direct summation of the Fresnel-Kirchoff
     * diffraction integral:
     *
     * \f{equation}{
     * U(x_1,y_1,z) = \frac{k}{2\pi i z} \int\int U(x,y,0) \exp\left\{ i k
     * \frac{ \left(x - x_1\right)^2 + \left(y - y_1\right)^2 } { 2z } \right\}
     * dx dy
     * \f}
     *
     * with functions \f$U(x,y,0)\f$ and \f$U(x,y,z)\f$ defined on rectangular grids. This
     * integral may be converted into a convolution form which can be
     * efficiently computed using FFT [7, 8]. This method is free from many
     * drawbacks of the spectral method given by the sequence
     * \f$(3) \rightarrow (2) \rightarrow (4)\f$, although it
     * is still very fast due to its use of FFT for computing of the integral
     * sums.
     *
     * We'll explain this using two-dimensional example, following [7], p.100.
     * Let the integral is defined in a finite interval
     * \f$-L/2 \ldots L/2\f$:
     *
     * \f{equation}{
     * U(x_1,z) = \sqrt{\frac{k}{2\pi i z}} \int_{-L/2}^{L/2} U(x,0) \exp
     * \left\{ i k \frac{ (x-x_1)^2 } {2z} \right\}dx
     * \f}
     *
     * Replacing functions \f$U(x)\f$ and \f$U(x_1)\f$ with step
     * functions \f$U_j\f$ and \f$U_m\f$,
     * defined in the sampling points of the grid with \f$j=0\ldots N\f$, and
     * \f$m=0\ldots N\f$ we convert the integral to the form:
     *
     * \f{equation}{
     * U_m = \sqrt{\frac{k}{2\pi i z}}\left( \sum_{j=1}^{N-1} U_j \int_{x_1 -
     * 0.5}^{x_1 + 0.5} \exp\left\{i k \frac{(x_m - x)^2}{2 z} \right\} dx + U_0
     * \int_{x_0}^{x_{0.5}} \exp\left\{i k \frac{(x_m - x)^2 }{2 z} \right\} dx
     * + U_N \int_{x_N - 0.5}^{X_N} \exp \left\{ i k \frac{ (x_m - x)^2}{2
     * z}\right\} dx \right)
     * \f}
     *
     * Performing these integrals, we obtain:
     *
     * \f{equation}{
     * U_m = \sum_{j=1}^{N-1} U_j K_{mj} + U_0 K_{m0} + U_N K_{mN}
     * \f}
     *
     * where: \f$K_{m0}, K_{mj}, K_{mN}\f$ are analytically expressed with the help
     * of Fresnel integrals, depending only onto the difference of indices. Sums
     * \f$\sum_{j=1}^{N-1} U_j K_{mj}\f$ can be easily calculated for all
     * indices \f$m\f$ as one convolution with the help of FFT.

     * The Fresnel filter implements this algorithm for two-dimensional
     * diffraction integrals. It is almost as fast as @ref forvard (still 2 to 5
     * times slower), it uses <b>8 times more memory</b> (to perform the
     * two-dimensional convolution) than @ref forvard and it allows for ``more
     * honest'' calculation of Fresnel diffraction.  As it does not require any
     * protection bands at the edges of the region, the model may be built in a
     * smaller grid, therefore the resources consumed and time of execution are
     * comparable or even better than that of @ref forvard. Fresnel does not
     * accept negative propagation distance. When possible Fresnel should be
     * used as the main propagation engine within LightPipes.
     *
     * Warning: @ref fresnel does not produce valid results if the distance of
     * propagation is comparable with (or less than) the characteristic size of
     * the aperture, at which the field is diffracted. In this case @ref forvard
     * or @ref steps should be used.
     *
     * @param z
     *     Distance to propagate.
     *
     * @see forvard, steps
     */
    Field & fresnel ( const double &z );

    /** Propagates Field using direct integration.
     *
     * Direct calculation of the Fresnel-Kirchoff integrals is very inefficient
     * in two-dimensional grids. The number of operations is proportional to
     * \f$N^4\f$, where \f$N\f$ is the grid sampling. With direct integration we
     * do not have any reflection at the grid boundary, so the size of the grid
     * can just match the cross section of the field distribution.  Forward has
     * following features:
     *    - arbitrary sampling and size of square grid at the input plane;
     *    - arbitrary sampling and size of square grid at the output plane, it
     *      means we can propagate field from a grid containing for example 52x52
     *      points corresponding to 4.9x4.9cm to a grid containing 42x42 points
     *      and corresponding let's say 8.75x8.75 cm.
     *
     * The direct calculation of diffraction integral in the Fresnel
     * approximation is used, thus forward can not be used to propagate the
     * field to a short distances.  Use forvard or steps instead.
     *
     * forward is a <b><i>very</i></b> slow routine.
     * @see forvard, fresnel, steps
     */
    Field & forward( const double & z,
                     const Pair<double> & new_side_length,
                     const Pair<size_t> & new_number );

    /** Propagates Field in spherical coordinates using FFT.
     * @param f
     *     Focal distance of lens that determines the curvature of the
     *     coordinate system. 
     * @param z
     *     Propagation distance. 
     */
    Field & lens_forvard(double f, double z);

    /** Propagates Field in spherical coordinates.
     * @param f
     *     Focal distance of lens that determines the curvature of the
     *     coordinate system. 
     * @param z
     *     Propagation distance. 
     */
    Field & lens_fresnel ( const double & f, const double & z );

    /** Propagates Field using FFT. 
     * @param z
     *     Distance to propagate.
     */
    Field & forvard( const double & z );

    /** Convert from spherical coordinates to normal coordinates. */
    Field & spherical_to_normal_coords ( );

    Field & circular_aperture( const double & r,
                               const double & x0 = 0.0,
                               const double & y0 = 0.0 );

    Field & circular_screen( const double & r,
                             const double & x0 = 0.0,
                             const double & y0 = 0.0 );

    /** Cut a rectangular aperture in the field.  zero all field outside of
     * aperture.
     * @param Lx  Size in x-direction
     * @param Ly  Size in y-direction.  If Ly < 0, use Ly = Lx
     * @param x0  x-coordinate of center of rectangle.
     * @param y0  y-coordinate of center of rectangle.
     * @param angle  angle of axes of rectangle.
     * @return    Reference to this field structure.
     */
    Field & rectangular_aperture ( const double & Lx,
                                         double   Ly = -0.1,
                                   const double & x0 = 0.0,
                                   const double & y0 = 0.0,
                                   const double & angle = 0.0 ); 

    /** Mask a rectangular aperture in the field.  zero all field outside of
     * aperture.
     * @param Lx  Size in x-direction
     * @param Ly  Size in y-direction.  If Ly < 0, use Ly = Lx
     * @param x0  x-coordinate of center of rectangle.
     * @param y0  y-coordinate of center of rectangle.
     * @param angle  angle of axes of rectangle.
     * @return    Reference to this field structure.
     */
    Field & rectangular_screen   ( const double & Lx,
                                         double   Ly = -0.1,
                                   const double & x0 = 0.0,
                                   const double & y0 = 0.0,
                                   const double & angle = 0.0 ); 

    /** Supergaussian aperature.
     * @param w
     *    1/e intensity width
     * @param n
     *    (2*n) == power of super gaussian
     */
    Field & supergaussian_aperture(const double & w,
                                   const int    & n,
                                   const double & x0 = 0.0,
                                   const double & y0 = 0.0,
                                   const double & A = 1.0 );

    /** Supergaussian aperature.
     * @param w
     *    1/e intensity width
     * @param n
     *    (2*n) == power of super gaussian
     */
    Field & supergaussian_screen  (const double & w,
                                   const int    & n,
                                   const double & x0 = 0.0,
                                   const double & y0 = 0.0,
                                   const double & A = 1.0 );

    Field & gaussian_aperture    ( const double & w,
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & A = 1.0 ) {
      return supergaussian_aperture( w, 1, x0, y0, A );
    }

    Field & gaussian_screen      ( const double & w,
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & A = 1.0 ) {
      return supergaussian_screen( w, 1, x0, y0, A );
    }

    Field & fft3 ( int ind );

    Field & tilt ( double tx, double ty );

    /** Zernike polynomials.
     * @param n
     *  Radial order
     * @param m
     *  Azimuthal order
     * @param R
     *  Radius of amplitute spec
     * @param A
     *  Amplitude of abberation at r=R, in radians
     */
    Field & zernike ( int n, int m, double R, double A );

    std::ostream & print_strehl (std::ostream & output) const ;

    std::pair<double,double> get_strehl_and_energy () const ;

    double get_strehl () const ;

    Field& pip_fft(const int&);

    /** Test to see if two Fields are compatible based on their field type
     * (wavelength, size, etc.).
     */
    bool compatible (const Field & that) const { return info.compatible(that.info); }

    /** Normalize the field to have unity total power. 
     * @param norm_coeff
     *    Output variable of the normalization factor applied. 
     */
    Field & normalize(double * norm_coeff = NULL );
    Field & l_amplify(const double & gain, const double & length, const double & i_sat);

    /** Interpolate the field onto a new grid.
     * @param angle Angle of rotation of field in radians. 
     */
    Field & interpolate(const Pair<double> & new_side_length = 0.0,
                        const Pair<size_t> & new_number = 0,
                        const double & x_shift = 0.0,
                        const double & y_shift = 0.0,
                        const double & angle = 0.0,
                        const double & magnif = 1.0) throw (std::runtime_error);


    /** Propagate the field using finite-difference routine.
     *
     */
    Field & steps( const double & step_size,
                   const int & number_steps,
                   const std::string & n_filename,
                   const std::string & k_filename = "",
                   const std::string & dump_filename = "",
                   const int & dump_period = 1 ) throw (std::runtime_error);

    /** Propagate the field using finite-difference routine.
     * @param n  Must be of size info.number+2(?) if not NULL
     */
    Field & steps( const double & step_size,
                   const int & number_steps = 1,
                   const Pixel * n = NULL,
                   const std::string & dump_filename = "",
                   const int & dump_period = 1 ) throw (std::runtime_error);

    /* **** END FIELD PHYSICAL OPERATORS. ***** */



    const Field::Info & getinfo() const { return info; }




    /** Print norm values of Field.  Calls print_field. */
    std::ostream & print_norm(std::ostream & output,
                               Pair<size_t> output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false) {
      return print_field(output,output_size,gamma,max_val,ascii,true);
    }

    /** Print phase values of Field. Calls print_field. */
    std::ostream & print_phase(std::ostream & output,
                               Pair<size_t> output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false) {
      return print_field(output,output_size,gamma,max_val,ascii,false);
    }

    /** Print either norm or phase values of Field. */
    std::ostream & print_field(std::ostream & output,
                               Pair<size_t> output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false,
                               const bool & output_norm = true);



  private:
    void cleanup();
    void base_init();
    void init(Pixel init_fill);
  };


  /** Field addition operator. */
  inline Field operator + (const Field & f1, const Field & f2) throw (std::runtime_error)  {
    Field retval(f1);
    retval += f2;
    return retval;
  }

}/* namespace lightpipes */

#endif // lightpipies_Field_h
