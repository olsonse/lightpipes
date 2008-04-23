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

#ifndef FIELD_H
#define FIELD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdexcept>
#include <complex>

#ifdef _DJGPP_
#  include <fcntl.h>
#endif

#ifndef SQR
template <class T> T SQR(const T & t) { return t*t; }
#endif



#ifdef _DJGPP_
extern void setmode ( int, int );
#endif


/*
 * The structure FIELD contains the characteristics of the light beam:
 * number of points along the side of a square grid, wavelength and side
 * length of the square grid, then two huge arrays of Re and Im data 
 */

class Field {
  private:
    void cleanup() { if (val) { delete[] val; val = NULL; } }

    void init () {
        cleanup();
        val = new std::complex<double>[SQR(info.number)];
        if (val == NULL)
            throw std::runtime_error("field allocation error");
        memset(val, 0, sizeof(std::complex<double>) * SQR(info.number));
    }

  public:
    class Info {
      public:
        Info() : number(0), side_length(0.), lambda(0.),
                 fft_level(0), sph_coords_factor(0.) { }
        int number;
        double side_length, lambda;
        int fft_level;
        double sph_coords_factor;

        bool compatible(const Info & that) const {
            if (number  != that.number  ||
                side_length    != that.side_length    ||
                lambda  != that.lambda  ||
                sph_coords_factor != that.sph_coords_factor)
                return false;
            else return true;
        }

        std::istream & read(std::istream & in) {
            if ( !in.read((char*)&number, sizeof(int)) )
                throw std::runtime_error("Error while reading FIELD.number\n");

            if ( number / 2 != ( float ) ( number ) / 2. )
                throw std::runtime_error("Sorry, number of points must be even, stopping\n");


            if ( !in.read((char*)&side_length, sizeof(double)) )
                throw std::runtime_error("Error while reading FIELD.side_length\n");

            if ( !in.read((char*)&lambda, sizeof(double)) )
                throw std::runtime_error("Error while reading FIELD.lambda\n");

            if ( !in.read((char*)&fft_level, sizeof(int)) )
                throw std::runtime_error("Error while reading FIELD.fft_level\n");

            if ( !in.read((char*)&sph_coords_factor, sizeof(double)) )
                throw std::runtime_error("Error while reading FIELD.sph_coords_factor\n");

            return in;
        }

        std::ostream & write(std::ostream & out) {
            if ( !out.write((char*)&number, sizeof(int)) )
                throw std::runtime_error("Error while writing FIELD.number\n");

            if ( !out.write((char*)&side_length, sizeof(double)) )
                throw std::runtime_error("Error while writing FIELD.side_length\n");

            if ( !out.write((char*)&lambda, sizeof(double)) )
                throw std::runtime_error("Error while writing FIELD.lambda\n");

            if ( !out.write((char*)&fft_level, sizeof(int)) )
                throw std::runtime_error("Error while writing FIELD.fft_level\n");

            if ( !out.write((char*)&sph_coords_factor, sizeof(double)) )
                throw std::runtime_error("Error while writing FIELD.sph_coords_factor\n");

            return out;
        }
    };


    /** Constructor. */
    Field ( unsigned int number, double side_length, double lambda,
            int fft_level = 0, double sph_coords_factor = 0.0 ) {
        val = NULL;
        info.number = number;
        info.side_length = side_length;
        info.lambda = lambda;
        info.fft_level = fft_level;
        info.sph_coords_factor = sph_coords_factor;
        init();
    }

    Field (const Info & that_info) {
        val = NULL;
        info = that_info;
        init();
    }

    Field (const Field & that) {
        val = NULL;
        (*this) = that;
    }

    ~Field () { cleanup(); }

    std::ostream & write(std::ostream & out = std::cout);
    static Field * read(std::istream & in = std::cin) throw (std::runtime_error) ;

    std::complex<double> & operator[](const unsigned int i) { return val[i]; }
    const std::complex<double> & operator[](const unsigned int i) const { return val[i]; }

    Field & operator = (const Field & that) {
        cleanup();
        this->info = that.info;

        val = new std::complex<double>[SQR(info.number)];
        if (val == NULL)
            throw std::runtime_error("field allocation error");
        memcpy(val, that.val, sizeof(std::complex<double>) * SQR(info.number));
        return *this;
    }

    template <class T>
    Field & operator *= ( const T & input ) {
        for (int i=1;i<=info.number; i++){
            for (int j=1;j<=info.number; j++){
                long ik1=(i-1)*info.number+j-1; 
                val[ik1] *= input;
            }
        }

        return *this;
    }

    Field & operator += ( const Field & that ) {
        if (!compatible(that))
            throw std::runtime_error("cannot add fields that are not compatible");

        for ( int i = SQR(info.number)-1; i >= 0; i-- ) {
            val[i] += that.val[i];
        }

        return *this;
    }


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
    Field & t_lens ( const double & R, const double & f, const double & x0, const double & y0 );

    /** Apply an axicon lens.
     * @param phi including angle of axicon (note that the sign of this value
     * is disregarded).
     * @param n1 complex index of refraction of axicon material.
     * @param x0 x-shift in position of center.
     * @param y0 y-shift in position of center.
     */
    Field & axicon ( const double & phi, const std::complex<double> & n1, const double & x0, const double & y0 );

    Field & fresnel ( const double &z );

    /**
     * Note that this operates on the input field structure.
     * @returns a reference to the changed field.
     */
    Field & forward( const double & z, const double & new_side_length, const int & new_number );
    Field & lens_forvard(double f, double z);
    Field & lens_fresnel ( const double & f, const double & z );
    Field & forvard( const double & z );
    Field & spherical_to_normal_coords ( );

    Field & circular_aperture( const double & r, const double & x0, const double & y0);
    Field & circular_screen( const double & r, const double & x0 = 0.0, const double & y0 = 0.0 );
    Field & rectangular_aperture ( const double & Lx, const double & Ly = -0.1,
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & angle = 0.0 ); 
    Field & rectangular_screen   ( const double & Lx, const double & Ly = -0.1,
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & angle = 0.0 ); 

    Field & supergaussian_aperture(const double & w, /* 1/e intensity width */
                                   const int    & n, /* (2*n) == power of super gaussian */
                                   const double & x0 = 0.0, const double & y0 = 0.0,
                                   const double & A = 1.0 );
    Field & supergaussian_screen  (const double & w, /* 1/e intensity width */
                                   const int    & n, /* (2*n) == power of super gaussian */
                                   const double & x0 = 0.0, const double & y0 = 0.0,
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

    bool compatible (const Field & that) const { return info.compatible(that.info); }

    Field & normalize(double * norm_coeff = NULL );
    Field & l_amplify(const double & gain, const double & length, const double & i_sat);

    std::ostream & print_field(std::ostream & output,
                               int output_size = 0,
                               const double & gamma = 2.0,
                               const int & max_val = 255,
                               const bool ascii = false);

    Info info;
    std::complex < double > * val;
};


inline Field operator + (const Field & f1, const Field & f2) throw (std::runtime_error)  {
    Field retval(f1);
    retval += f2;
    return retval;
}

#endif // FIELD_H
