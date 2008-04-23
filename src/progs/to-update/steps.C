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
#include <math.h>
#include <lightpipes/Field.h>

void all_mem (  ), free_mem (  ), elim (  );

void error_print ( const char * arr = NULL );

std::complex<double> *a, *b, *c, *alpha, *beta, *u, *p, *u1, *u2;
double *refr, *absorb;

int main ( int argc, char *argv[] ) {
    void step (  );
    long ik, ik1, ikij, ikij1, ikij_1, iki1j, iki_1j;
    double z, delta, delta2, Pi4lz, A, AA, band_pow, K, dx, dist;
    std::complex<double> uij, uij1, uij_1, ui1j, ui_1j, medium;
    int i, j, jj, ii, nstep, istep, i_left, i_right, n_st;
    FILE *fr = NULL;
    double *int1 = NULL;
    double int2;
    double *phase1 = NULL;
    double phase2;

    /*
     * Processing the command line argument 
     */
    if ( argc < 2 || argc > 7 ) {
        error_print ( argv[0] );
        exit ( 1 );
    }

    Field * field = Field::read (  );

    all_mem (  );

    K = 2. * Pi / info.lambda;
    /*
     * reading the data from a command line 
     */
    sscanf ( argv[1], "%le", &z );
    nstep = 1;
    if ( argc > 2 ) sscanf ( argv[2], "%d", &nstep );

    if ( argc > 3 ) {
        if ( ( strstr ( argv[3], "void" ) ) != NULL ) {
            fprintf ( stderr,
                      "%s: void refractive coefficient file, skipping \n",
                      argv[0] );
            goto l1;
        }

        if ( ( fr = fopen ( argv[3], "r" ) ) == NULL ) {
            fprintf ( stderr, "%s: error opening file %s, exiting \n",
                      argv[0], argv[3] );
            exit ( 1 );
        }

        ik = 0;
        for ( i = 1; i <= info.number; i++ ) {
            for ( j = 1; j <= info.number; j++ ) {
                double fi;
                if ( ( fscanf ( fr, "%le", &fi ) ) == EOF ) {
                    fprintf ( stderr, "%s: reading the refractive indices:\
  end of input file reached, exiting\n", argv[0] );
                    exit ( 1 );
                }

                refr[ik] = fi;
                ik++;

            }
        }
        fclose ( fr );
    }

  l1:
    if ( argc > 4 ) {
        if ( ( strstr ( argv[4], "void" ) ) != NULL ) {
            fprintf ( stderr,
                      "%s: void absorption coefficient file, skipping \n",
                      argv[0] );
            goto l2;
        }

        if ( ( fr = fopen ( argv[4], "r" ) ) == NULL ) {
            fprintf ( stderr, "%s: error opening file %s, exiting \n",
                      argv[0], argv[4] );
            exit ( 1 );
        }

        ik = 0;
        for ( i = 1; i <= info.number; i += 1 ) {
            for ( j = 1; j <= info.number; j += 1 ) {
                double fi;
                if ( ( fscanf ( fr, "%le", &fi ) ) == EOF ) {
                    fprintf ( stderr,
                              "%s: absorption  end of input file reached, exiting\n",
                              argv[0] );
                    exit ( 1 );
                }
                absorb[ik] = fi;
                ik++;

            }
        }

        fclose ( fr );
    }

  l2:

    if ( argc > 5 ) {

#ifdef _DJGPP_
        setmode ( fileno ( stdout ), O_BINARY );

        fr = fopen ( argv[5], "wb" );
#else
        fr = fopen ( argv[5], "w" );
#endif

        int1 = new double[info.number + 2];
        memset ( int1, 0, info.number + 2 );

        if ( int1 == NULL ) {
            fprintf ( stderr, "Allocation error in %s , int1\n", argv[0] );
            exit ( 1 );
        }

        phase1 =
            ( double * ) calloc ( info.number + 2, sizeof ( double ) );
        if ( phase1 == NULL ) {
            fprintf ( stderr, "Allocation error in %s , phase1\n",
                      argv[0] );
            exit ( 1 );
        }

    }
    dx = info.side_length / ( info.number - 1. );
    n_st = 1;

    if ( argc > 6 ) {
        sscanf ( argv[6], "%d", &n_st );
    }


    /*
     * the arrays for refraction and absorption are finally formed 
     */




    if ( info.sph_coords_factor != 0. ) {
        fprintf ( stderr, "%s can not be applied in spherical"
                          " coordinates,\nuse CONVERT first\n", argv[0] );
        exit ( 1 );
    }


    z = z / 2.;
    Pi4lz = 4. * Pi / info.lambda / z;
    delta = info.side_length / ( ( double ) ( info.number - 1. ) );
    delta2 = delta * delta;
    A = 0.;
    /*
     * absorption at the borders is described here 
     */
    AA = -10. / z / nstep;      /* total absorption */
    band_pow = 2.;              /* profile of the absorption border,
                                 * 2=quadratic */
    /*
     * width of the absorption border 
     */
    i_left = info.number / 2 + 1 - 0.4 * info.number;
    i_right = info.number / 2 + 1 + 0.4 * info.number;
    /*
     * end absorption 
     */


    for ( i = 1; i <= info.number; i++ ) {
        u2[i].r = u2[i].i = 0.;
        /*
         * a[i].r=b[i].r= -1./delta2; 
         */
        a[i].i = b[i].i = 0.;
        a[i].r = b[i].r = -1. / delta2;
    }

    medium.r = 0.;
    dist = 0.;
    /*
     * Main loop, steps here 
     */
    for ( istep = 1; istep <= nstep; istep++ ) {
        dist = dist + 2. * z;
        /*
         * Elimination in the direction i, halfstep 
         */

        ik = 0;
        for ( i = 1; i <= info.number; i++ ) {
            for ( j = 1; j <= info.number; j++ ) {
                register double fi = 0.25 * K * z * ( refr[ik] - 1. );
                val[ik] *= exp(I * fi);

                ik++;
            }
        }


        for ( jj = 2; jj <= info.number - 2; jj += 2 ) {

            j = jj;
            for ( i = 2; i <= info.number - 1; i++ ) {

                ikij = ( j - 1 ) * info.number + i - 1;
                ikij1 = ( j ) * info.number + i - 1;
                ikij_1 = ( j - 2 ) * info.number + i - 1;


                uij = val[ikij];
                uij1 = val[ikij1];
                uij_1 = val[ikij_1];

                p[i] = -2. * uij;
                p[i] = p[i] + uij1;
                p[i] = p[i] + uij_1;
                p[i] = ( -1. / delta2 ) * p[i];
                p[i] =  std::complex<double> ( 0., Pi4lz ) * uij +  p[i];

            }

            for ( i = 1; i <= info.number; i++ ) {
                ik = ( j - 1 ) * info.number + i - 1;


                medium.i = -2. * Pi * absorb[ik] / info.lambda;
                c[i].r = -2. / delta2;
                c[i].i = Pi4lz;
                c[i].i += medium.i;
                /*
                 * absorption borders are formed here 
                 */
                if ( i <= i_left ) {
                    int iii = i_left - i;
                    c[i].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( i_left ) ),
                                           band_pow );
                }
                if ( i >= i_right ) {
                    int iii = i - i_right;
                    int im = info.number - i_right;
                    c[i].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( im ) ),
                                           band_pow );
                }
                /*
                 * end absorption 
                 */

            }
            elim (  );

            for ( i = 1; i <= info.number; i++ ) {
                ikij_1 = ( j - 2 ) * info.number + i - 1;
                val[ikij_1] = u2[i];
                u1[i] = u[i];
            }

            j = jj + 1;

            for ( i = 2; i <= info.number - 1; i++ ) {

                ikij = ( j - 1 ) * info.number + i - 1;
                ikij1 = ( j ) * info.number + i - 1;
                ikij_1 = ( j - 2 ) * info.number + i - 1;


                uij = val[ikij];
                uij1 = val[ikij1];
                uij_1 = val[ikij_1];

                p[i] = -2. * uij;
                p[i] = p[i] + uij1;
                p[i] = p[i] + uij_1;
                p[i] = ( -1. / delta2) * p[i];
                p[i] = std::complex<double> ( 0., Pi4lz ) * uij + p[i];

            }
            for ( i = 1; i <= info.number; i++ ) {
                ik = ( j - 1 ) * info.number + i - 1;


                medium.i = -2. * Pi * absorb[ik] / info.lambda;
                c[i].r = -2. / delta2;
                c[i].i = Pi4lz;
                c[i].i += medium.i;

                /*
                 * absorption borders are formed here 
                 */
                if ( i <= i_left ) {
                    int iii = i_left - i;
                    c[i].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( i_left ) ),
                                           band_pow );
                }
                if ( i >= i_right ) {
                    int iii = i - i_right;
                    int im = info.number - i_right;
                    c[i].i -=
                        ( AA * 2. * K ) * pow ( ( double ) iii /
                                                ( ( double ) ( im ) ),
                                                band_pow );
                }
                /*
                 * end absorption 
                 */

            }
            elim (  );

            for ( i = 1; i <= info.number; i++ ) {
                ikij = ( j - 2 ) * info.number + i - 1;
                val[ikij] = u1[i]
                u2[i] = u[i];
            }
        }

        for ( i = 1; i <= info.number; i++ ) {
            ikij = ( info.number - 1 ) * info.number + i - 1;
            val[ikij] = u2[i];

        }
        ik = 0;
        for ( i = 1; i <= info.number; i++ ) {
            for ( j = 1; j <= info.number; j++ ) {
                register double fi = 0.5 * K * z * ( refr[ik] - 1. );
                val[ik] *= exp(I * fi);

                ik++;
            }
        }

        /*
         * Elimination in the j direction is here, halfstep 
         */

        for ( i = 1; i <= info.number; i++ ) {
            u2[i].r = u2[i].i = 0.;
        }

        for ( ii = 2; ii <= info.number - 2; ii += 2 ) {

            i = ii;
            for ( j = 2; j <= info.number - 1; j++ ) {

                ikij = ( j - 1 ) * info.number + i - 1;
                iki1j = ( j - 1 ) * info.number + i;
                iki_1j = ( j - 1 ) * info.number + i - 2;


                uij = val[ikij];
                ui1j = val[iki1j];
                ui_1j = val[iki_1j];

                p[j] = -2. * uij;
                p[j] = p[j] + ui1j;
                p[j] = p[j] + ui_1j;
                p[j] = ( -1. / delta2 ) * p[j];
                p[j] = std::complex<double> ( 0., Pi4lz ) * uij + p[j];

            }
            for ( j = 1; j <= info.number; j++ ) {
                ik = ( j - 1 ) * info.number + i - 1;


                medium.i = -2. * Pi * absorb[ik] / info.lambda;
                c[j].r = -2. / delta2;
                c[j].i = Pi4lz;
                c[j].i += medium.i;

                /*
                 * absorption borders are formed here 
                 */
                if ( j <= i_left ) {
                    int iii = i_left - j;
                    c[j].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( i_left ) ),
                                           band_pow );
                }
                if ( j >= i_right ) {
                    int iii = j - i_right;
                    int im = info.number - i_right;
                    c[j].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( im ) ),
                                           band_pow );
                }
                /*
                 * end absorption 
                 */

            }
            elim (  );

            for ( j = 1; j <= info.number; j++ ) {
                iki_1j = ( j - 1 ) * info.number + i - 2;
                val[iki_1j] = u2[j];
                u1[j] = u[j];
            }

            i = ii + 1;
            for ( j = 2; j <= info.number - 1; j++ ) {

                ikij = ( j - 1 ) * info.number + i - 1;
                iki1j = ( j - 1 ) * info.number + i;
                iki_1j = ( j - 1 ) * info.number + i - 2;


                uij = val[ikij];
                ui1j = val[iki1j];
                ui_1j = val[iki_1j];

                p[j] = -2. * uij;
                p[j] = p[j] + ui1j;
                p[j] = p[j] + ui_1j;
                p[j] = ( -1. / delta2 ) * p[j];
                p[j] = std::complex<double> ( 0., Pi4lz ) * uij + p[j];

            }
            for ( j = 1; j <= info.number; j++ ) {
                ik = ( j - 1 ) * info.number + i - 1;


                medium.i = -2. * Pi * absorb[ik] / info.lambda;
                c[j].r = -2. / delta2;
                c[j].i = Pi4lz;
                c[j].i += medium.i;

                /*
                 * absorption borders are formed here 
                 */
                if ( j <= i_left ) {
                    int iii = i_left - j;
                    c[j].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( i_left ) ),
                                           band_pow );
                }
                if ( j >= i_right ) {
                    int iii = j - i_right;
                    int im = info.number - i_right;
                    c[j].i -=
                        ( AA * K ) * pow ( ( double ) iii /
                                           ( ( double ) ( im ) ),
                                           band_pow );
                }
                /*
                 * end absorption 
                 */

            }
            elim (  );


            for ( j = 1; j <= info.number; j++ ) {
                ikij = ( j - 1 ) * info.number + i - 2;
                val[ikij] = u1[j];
                u2[j] = u[j];
            }
        }

        for ( j = 1; j <= info.number; j++ ) {
            ikij = ( j - 1 ) * info.number + i - 1;
            val[ikij] = u2[j];
        }

        /*
         * end j 
         */

    /***************************************************/


        if ( ( argc > 5 )
             && ( istep / n_st == ( float ) istep / ( float ) n_st ) ) {

            /*
             * writing the intensity into arrays 
             */

            i = info.number / 2 + 1;
            for ( j = 1; j <= info.number; j += 1 ) {
                ik1 = ( i - 1 ) * info.number + j - 1;
                jj = j - 1;
                int1[jj] = norm( val[ik1] );
                phase1[jj] = arg(val[ik1]);
            }


            j = info.number / 2 + 1;
            for ( i = 1; i <= info.number; i += 1 ) {
                double cc;
                ik1 = ( i - 1 ) * info.number + j - 1;
                jj = i - 1;
                cc = dx * ( i - info.number / 2 - 1 );

                int2 = norm (val[ik1]);
                phase2 = arg (val[ik1]);
                fprintf ( fr, " %e %e %e %e %e %e\n", cc, int1[jj], int2,
                          phase1[jj], phase2, dist );

            }


        }
        if ( ( argc > 5 )
             && ( istep / n_st == ( float ) istep / ( float ) n_st ) ) {
            fprintf ( fr, "\n" );
        }

    }


    ik = 0;
    for ( i = 1; i <= info.number; i++ ) {
        for ( j = 1; j <= info.number; j++ ) {
            register double fi = 0.25 * K * z * ( refr[ik] - 1. );
            val[ik] *= exp(I * fi);

            ik++;
        }
    }

    free_mem (  );

    if ( argc > 5 ) {
        free ( int1 );
        free ( phase1 );
        fclose ( fr );
    }

    field->write (  );
    delete field;

    return 0;
}


void elim (  ) {
    int i;
    std::complex<double> cc;
    /*
     * initial condition, everything is going to be zero at the edge 
     */
    alpha[2].r = 0.;
    alpha[2].i = beta[2].r = beta[2].i = 0.;

    alpha[info.number].r = 0.;
    alpha[info.number].i = beta[info.number].r = beta[info.number].i =
        0.;

    /*
     * forward elimination 
     */
    for ( i = 2; i <= info.number - 2; i++ ) {
        cc = c[i] -  a[i] * alpha[i];
        alpha[i + 1] = b[i] / cc;
        beta[i + 1] = ( p[i] + a[i] * beta[i] ) / cc;
    }
    i = info.number;
    cc = c[i] - a[i] * alpha[i];
    beta[info.number + 1] = ( p[i] + ( a[i] * beta[i] ) ) / cc;
    /*
     * edge amplitude =0 
     */
    u[info.number] = beta[info.number + 1];
    /*
     * backward elimination 
     */
    for ( i = info.number - 1; i >= 1; i-- ) {
        u[i] = ( alpha[i + 1] * u[i + 1] ) + beta[i + 1];
    }
}

void all_mem (  ) {
    int i, j;
    long ik;
    /*
     * Allocatiom of the memory for arrays: 
     */

    if ( ( refr =
           ( double * ) calloc ( ( info.number + 2 ) *
                                 ( info.number + 2 ),
                                 sizeof ( double ) ) ) == NULL ) {
        fprintf ( stderr,
                  "Steps: allocation error, array of refractive indices" );
        exit ( 1 );
    }

    if ( ( absorb =
           ( double * ) calloc ( ( info.number + 2 ) *
                                 ( info.number + 2 ),
                                 sizeof ( double ) ) ) == NULL ) {
        fprintf ( stderr, "Steps allocation error, array of absorption" );
        exit ( 1 );
    }

    ik = 0;
    for ( i = 1; i <= info.number + 2; i++ ) {
        for ( j = 1; j <= info.number + 2; j++ ) {
            refr[ik] = 1.;
            absorb[ik] = 0.;
            ik++;
        }
    }

    /*
     * end arrays 
     */
    /*
     * allocation of memory for the elimination 
     */

    if ( ( a =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    if ( ( p =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    if ( ( b =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }
    if ( ( c =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    if ( ( u =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    if ( ( u1 =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    if ( ( u2 =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }
    if ( ( alpha =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }
    if ( ( beta =
           ( std::complex<double> * ) calloc ( info.number + 3,
                                   sizeof ( std::complex<double> ) ) ) == NULL ) {
        fprintf ( stderr, "Allocation error" );
        exit ( 1 );
    }

    /*
     * Memory done 
     */

}

void free_mem (  ) {
    /*
     * freeing the memory 
     */

    free ( a );
    free ( b );
    free ( c );
    free ( u );
    free ( u1 );
    free ( u2 );
    free ( alpha );
    free ( beta );
    free ( p );
    free ( absorb );
    free ( refr );

}


/*
 * complex number handling is here, not optimal ... : 
 */

void error_print ( const char *arr ) {

    fprintf ( stderr, "\n%s propagates the field to\n\
distance N*Z [units you use]\n\
using finite difference method\n", arr );

    fprintf ( stderr, "\nUSAGE: \n" );

    fprintf ( stderr, "%s Z [N R A F M] ,\n\
where Z is the step size, N is the number of steps to perform\n\
R is the name of the file, which  contains the distribution\n\
of refractive index in gnuplot format\n\
if R equals to < void > then this option is skipped\n\
A is the name of a file which  contains the two-dimensional\n\
distribution of the absorption coefficient\n\
if A equals to < void > then this is skipped\n\
F is the filename where cross section of the beam\n\
is written at each M-th step\n\n", arr );

    fprintf ( stderr, "Examples:\n\
%s 0.1 20 - performs 20 steps of 0.1 m in a free space,\n\
total distance of propagation is  2 m \n\n\
%s 0.1 20 void abs_d - performs 20 steps of 0.1 m in a medium\n\
absorption coefficient of which is given in a file abs_d,\n\
refraction coefficient equals to unity\n\
%s 0.1 20 void void file_out - performs 20 steps\n\
in a free space writing cross section of the beam\n\
into the file file_out, in the format:\n\n\
x  I(x,0) I(0,y ) F(x,0) F(0,y) Z \n\n\
where I is the intensity and F is phase-\n\
the same format as used by cros_out.\n\
Three-dimensional profile of propagating beam\n\
can be plotted for example with gnuplot commands:\n\
set par; splot 'file_out' using 1:6:2 w l;\n\n", arr, arr, arr );


    fprintf ( stderr, "%s 0.1 200 refr_d abs_d out.dat 10\n\
- performs 200 steps of 0.1 m in a medium\n\
absorption coefficient of which is given in a file abs_d,\n\
refraction coefficient is given in a file refr_d\n\
the profile of the beam is written into file out.dat at each 10-th step\n\n\
%s 0.1 20 refr_d - performs 20 steps of 0.1 m in a medium\n\
refraction coefficient of which is given in a file refr_d\n\
without absorption\n\n", arr, arr );

}
