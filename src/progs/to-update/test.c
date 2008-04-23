/* check that exp(i*pi) == -1 */
#include <math.h>   /* for atan */
#include <complex.h>
main() {
    double pi = 4*atan(1);
    complex z = cexp(I*pi);
    printf("%f+%f*i\n", creal(z), cimag(z));
}
