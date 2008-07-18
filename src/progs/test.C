/* check that exp(i*pi) == -1 */
#include <cmath>   /* for atan */
#include <complex>
#include <iostream>
#include <stdexcept>

/** Gleb Vdovin's way of computing phase. */
double phase ( double y, double x ) throw (std::runtime_error) {
    double pp = 0.;
    if ( x == 0. ) {
        if ( y > 0  ) pp = 0.5 * M_PI;
        if ( y == 0 ) pp = 0.;
        if ( y < 0  ) pp = -0.5 * M_PI;
    } else {
        if ( y != 0 ) pp = atan2 ( y, x );
        else          pp = 0.;
    }
    return pp;
}


/** program to compare Vdovin's phase() function with
 * std::arg(std::complex<> &) function. 
 */
int main() {
    double pi = 4*atan(1);
    std::complex<double> I(0,1);
    std::complex<double> z = 2.*exp(I*pi/4.0);
    std::cout << "phase () : " << phase(z.imag(), z.real()) << std::endl
              << "arg ()   : " << arg(z) << std::endl
              << "abs ()   : " << abs(z) << std::endl
              << "norm()   : " << norm(z) << std::endl
              << "conj()   : " << conj(z) << std::endl
              << "z        : " << z << std::endl;

    int i = 0;
    for (double phi = 0; phi <= 10*M_PI; phi += 1e-6) {
        z.real() = cos(phi);
        z.imag() = sin(phi);
        double val_phase = phase(z.imag(), z.real());
        double val_arg = arg(z);
        if ( fabs(val_phase - val_arg)/abs(z) > 1e-6 ) {
            std::cout << "phase is different at " << z << std::endl;
            exit(EXIT_FAILURE);
        }
        
        i++;

        if ( fabs( z.imag() ) < 5e-7 ) {
            std::cout << "." << std::flush;
        }
    }

    std::cout << "tested " << i << " phase values " << std::endl;
}
