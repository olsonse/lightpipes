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
    /** Defines character of Field such as wavelength of light, spatial size,
     * spatial precision, etc. */
    struct Info {
      /* STORAGE MEMBERS */
      /** The number of side-elements stored in Field::val
       * (sizeof(Field::val) = number^{2}.
       */
      int number;

      /** The physical size of the side of the Field view. */
      double side_length;

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
    std::complex < double > * val;



    /* FUNCTION MEMBERS */
  public:
    /** Constructor. */
    Field ( unsigned int number,
            double side_length,
            double lambda,
            int fft_level = 0,
            double sph_coords_factor = 0.0 );

    /** Constructor with Field::Info initialization. */
    Field (const Info & that_info);

    /** Copy constructor. */
    Field (const Field & that);

    /** Destructor. */
    ~Field ();

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
    std::complex<double> & operator[](const unsigned int i) { return val[i]; }

    /** Field index operator (const). */
    const std::complex<double> & operator[](const unsigned int i) const { return val[i]; }

    /** Field index operator (non-const) using column AND row indices. */
    std::complex<double> & operator()(const unsigned int row,
                                      const unsigned int col) {
      return val[row*info.number + col];
    }

    /** Field index operator (const) using column AND row indices. */
    const std::complex<double> & operator()(const unsigned int row,
                                            const unsigned int col) const {
      return val[row*info.number + col];
    }

    /** Copy operator. */
    Field & operator= (const Field & that);

    /** fill the field with a constant value operator. */
    template < typename T >
    Field & fill(const T & t) {
      for (int i=1;i<=info.number; ++i){
        for (int j=1;j<=info.number; ++j){
          long ik1=(i-1)*info.number+j-1;
          val[ik1] = t;
        }
      }

      return *this;
    }

    /** Scaling operator. */
    template <class T>
    Field & operator*= ( const T & input ) {
      for (int i=1;i<=info.number; i++){
        for (int j=1;j<=info.number; j++){
          long ik1=(i-1)*info.number+j-1; 
          val[ik1] *= input;
        }
      }

      return *this;
    }

    /** Field addition operator (Field1 += Field2). */
    Field & operator+= ( const Field & that );


    /** Apply a lens.
     * f>0 for positive lens !!!!!  (Unlike the original lightpipes code)
     */
    Field & lens ( const double & f, const double & x0, const double & y0 );

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
                     const double & x0,
                     const double & y0 );

    /** Apply an axicon lens.
     * @param phi including angle of axicon (note that the sign of this value
     * is disregarded).
     * @param n1 complex index of refraction of axicon material.
     * @param x0 x-shift in position of center.
     * @param y0 y-shift in position of center.
     */
    Field & axicon ( const double & phi,
                     const std::complex<double> & n1,
                     const double & x0,
                     const double & y0 );

    /** Propagates Field using convolution.
     * @param z
     *     Distance to propagate.
     */
    Field & fresnel ( const double &z );

    /** Propagates Field using direct integration.
     * Note that this operates on the input field structure.
     * @returns a reference to the changed field.
     */
    Field & forward( const double & z,
                     const double & new_side_length,
                     const int & new_number );

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

    Field & rectangular_aperture ( const double & Lx,
                                   const double & Ly = -0.1,
                                   const double & x0 = 0.0,
                                   const double & y0 = 0.0,
                                   const double & angle = 0.0 ); 

    Field & rectangular_screen   ( const double & Lx,
                                   const double & Ly = -0.1,
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
      return supergaussian_aperture( w, 2, x0, y0, A );
    }

    Field & gaussian_screen      ( const double & w,
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & A = 1.0 ) {
      return supergaussian_screen( w, 2, x0, y0, A );
    }

    Field & fft3 ( int ind );

    Field & tilt ( double tx, double ty );

    Field & zernike ( int n, int m, double R, double A );

    std::ostream & print_strehl (std::ostream & output);

    double get_strehl ();

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
    Field & interpolate(const double & new_side_length = 0.0,
                        const int    & new_number = 0.0,
                        const double & x_shift = 0.0,
                        const double & y_shift = 0.0,
                        const double & angle = 0.0,
                        const double & magnif = 1.0) throw (std::runtime_error);


    /** Propagate the field using finite-difference routine.
     *
     */
    Field & steps( const double & step_size,
                   const int & N = 1,
                   const std::string & n_filename = "",
                   const std::string & k_filename = "",
                   const std::string & X_filename = "",
                   const int & dump_period = 1 ) throw (std::runtime_error);

    /* **** END FIELD PHYSICAL OPERATORS. ***** */



    const Field::Info & getinfo() const { return info; }




    /** Print norm values of Field.  Calls print_field. */
    std::ostream & print_norm(std::ostream & output,
                               int output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false) {
      return print_field(output,output_size,gamma,max_val,ascii,true);
    }

    /** Print phase values of Field. Calls print_field. */
    std::ostream & print_phase(std::ostream & output,
                               int output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false) {
      return print_field(output,output_size,gamma,max_val,ascii,false);
    }

    /** Print either norm or phase values of Field. */
    std::ostream & print_field(std::ostream & output,
                               int output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool & ascii = false,
                               const bool & output_norm = true);



  private:
    void cleanup();
    void init();
  };


  /** Field addition operator. */
  inline Field operator + (const Field & f1, const Field & f2) throw (std::runtime_error)  {
    Field retval(f1);
    retval += f2;
    return retval;
  }

}/* namespace lightpipes */

#endif // lightpipies_Field_h
