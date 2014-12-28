/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h>
#include "pipes.h"

void main(int argc, char *argv[]){

    double z, old_size, new_size, dx_new, dx_old, x_new, y_new,\
old_number, on21, nn21, *fieldimaginary, fieldreal, P1, P2, P3, P4, R22; 

double fc1, fs1, fc2, fs2, fc3, fs3, fc4, fs4, fr, fi;
double c4c1, c2c3, c4s1, s4c1, s2c3, c2s1, s4c3, s2c1, c4s3, s2s3, s2s1, c2s3, s4s1, c4c3, s4s3, c2c1;

    int new_number,i,j,i_old, j_old, io, jo, dum, fresnl();
    void error_print();
    long  kk, kk1;
    fs1=fc1=fs2=fc2=fs3=fc3=fs4=fc4=0.; /* to make the compiler happy */
    /* Processing the command line argument  */
#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);
#endif
    if (argc< 2 || argc >4 ){
        error_print(argv[0]);
        exit(1);
    }

    /* reading the data from a command line */

    sscanf(argv[1],"%le",&z);


    read_field();

    R22=sqrt(1./(2.*field.lambda*z));

    old_size=field.size;
    if (argc == 3) sscanf(argv[2],"%le",&new_size);
    else new_size=field.size;

    old_number=field.number;
    if (argc == 4){sscanf(argv[2],"%le",&new_size); 
     sscanf(argv[3],"%d",&new_number);}
    else new_number=field.number;
        
    dx_new=new_size/(new_number-1.); 
    dx_old=old_size/(old_number-1.);
    
   
    field.number=new_number;
    field.size=new_size;

    if(fwrite ((char *)&(field.number), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.number\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.size), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.size\n");
        exit(1);
    }
    if(fwrite ((char *)&(field.lambda), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.lambda\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int1), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.int1\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int2), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.int2\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int3), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.int3\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double1), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.double1\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double2), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.double2\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double3), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"forward: Error while writing FIELD.double3\n");
        exit(1);
    }

    on21=old_number/2+1;
    nn21=new_number/2+1;
 

/*  fieldimaginary is imaginary array in plane z */
    fieldimaginary = (double *) calloc( (new_number)*(new_number), sizeof(double) );
    if(fieldimaginary == NULL){fprintf(stderr,"\nAllocation error in forward\n");
			exit(1);}


/*           main cycle here starts            */

    kk=0;
   
    
    for (i=1; i<=field.number; i++){
        x_new=(i-nn21)*dx_new;
     
        for (j=1; j<=field.number; j++){
            y_new=(j-nn21)*dx_new;

fieldimaginary[kk]=0.;
fieldreal = 0.;

/* Inner cycle here - integral itself */

	    kk1=0;
  for (i_old=1; i_old<=old_number; i_old++){
        io=i_old-on21;
       
        for (j_old=1; j_old<=old_number; j_old++){
	  jo=j_old-on21;
          
	  P1=R22*(2*(dx_old*io-x_new)+dx_old);
	  P2=R22*(2*(dx_old*jo-y_new)-dx_old);
	  P3=R22*(2*(dx_old*io-x_new)-dx_old);
	  P4=R22*(2*(dx_old*jo-y_new)+dx_old);
	  dum=fresnl(P1,&fs1, &fc1);
	  dum=fresnl(P2,&fs2, &fc2);
	  dum=fresnl(P3,&fs3, &fc3);
	  dum=fresnl(P4,&fs4, &fc4);
	  fr=0.5*field.real[kk1];
	  fi=0.5*field.imaginary[kk1];
          c4c1=fc4*fc1;
	  c2s3=fc2*fs3;
	  c4s1=fc4*fs1;
	  s4c1=fs4*fc1;
	  s2c3=fs2*fc3;
	  c2s1=fc2*fs1;
	  s4c3=fs4*fc3;
	  s2c1=fs2*fc1;
	  c4s3=fc4*fs3;
fieldreal += fr*(c2s3+c4s1+s4c1+s2c3-c2s1-s4c3-s2c1-c4s3);
	  s2s3=fs2*fs3;
	  s2s1=fs2*fs1;
	  c2c3=fc2*fc3;
	  s4s1=fs4*fs1;
	  c4c3=fc4*fc3;
	  c4c1=fc4*fc1;
	  s4s3=fs4*fs3;
	  c2c1=fc2*fc1;
fieldreal += fi*(-s2s3+s2s1+c2c3-s4s1-c4c3+c4c1+s4s3-c2c1);

fieldimaginary[kk] += fr*(-c4c1+s2s3+c4c3-s4s3+c2c1-s2s1+s4s1-c2c3);
fieldimaginary[kk] += fi*(c2s3+s2c3+c4s1+s4c1-c4s3-s4c3-c2s1-s2c1); 
	    kk1 ++;
	  }
      }
          
	  

            if(fwrite (&fieldreal, sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"forward: error while writing FIELD.real  \n");
                exit (1);
            }
	    kk ++;
        }
    }

kk=0;
    for (i=1; i<=new_number; i++){
     
        for (j=1; j<=new_number; j++){
        
            if(fwrite (&fieldimaginary[kk], sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"forward: error while writing FIELD.imaginary  \n");
                exit (1);
            }
	    kk ++;
        }
    }

free(fieldimaginary);
}



void error_print(char *arr)
{
    fprintf(stderr,"\n%s propagates the field to \
distance Z [units you use]\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s  Z [x1 n1], where  Z is the distance to propagate\n\
x1 and n1 are the new grid size and number of points after propagation into plane Z\n\
direct calculation of the diffraction integral is implemented in forward\n\n",arr);

}

/* ----------Fresnel Integrals start here --------------- */
/* 
Fresnel Integrals from CEPHES  
Modified for portability by G. Vdovin (gleb@okotech.com     )
Oct 1995 
*/

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*/

/*							fresnl.c
 *
 *	Fresnel integrals
 *
 *
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

/*
This is a test program 

*/

/*

int fresnl();

main(){
double a,b,c;
int i, dum;
for (i=-1000; i<= 1000; i++){

a=(double)(((double) i)/100.);

dum=fresnl(a, &b, &c);

printf("%e %e %e \n", a,b,c);

}

}
*/



int fresnl( xxa, ssa, cca )
double xxa, *ssa, *cca;
{
double f, g, cc, ss, c, s, t, u;
double x, x2;
double polevl(), p1evl();

double PII= 3.14159265358979323844;
double PIO2 = 1.57079632679489661922;


static double sn[6] = {
-2.99181919401019853726E3,
 7.08840045257738576863E5,
-6.29741486205862506537E7,
 2.54890880573376359104E9,
-4.42979518059697779103E10,
 3.18016297876567817986E11,
};
static double sd[6] = {
/* 1.00000000000000000000E0,*/
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
/*  1.00000000000000000000E0,*/
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
/*  1.00000000000000000000E0,*/
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

x = fabs(xxa);
x2 = x * x;
if( x2 < 2.5625 )
	{
	t = x2 * x2;
	ss = x * x2 * polevl( t, sn, 5)/p1evl( t, sd, 6 );
	cc = x * polevl( t, cn, 5)/polevl(t, cd, 6 );
	goto done;
	}
if( x > 36974.0 )
	{
	cc = 0.5;
	ss = 0.5;
	goto done;
	}
/* Auxiliary functions for large argument  */
	x2 = x * x;
	t = PII * x2;
	u = 1.0/(t * t);
	t = 1.0/t;
	f = 1.0 - u * polevl( u, fn, 9)/p1evl(u, fd, 10);
	g = t * polevl( u, gn, 10)/p1evl(u, gd, 11);

	t = PIO2 * x2;
	c = cos(t);
	s = sin(t);
	t = PII * x;
	cc = 0.5  +  (f * s  -  g * c)/t;
	ss = 0.5  -  (f * c  +  g * s)/t;
done:
if( xxa < 0.0 )
	{
	cc = -cc;
	ss = -ss;
	}
*cca = cc;
*ssa = ss;
return(0);
}
/*							polevl.c
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


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


double polevl( x, coef, N )
double x;
double coef[];
int N;
{
double ans;
int i;
double *p;

p = coef;
ans = *p++;
i = N;

do
	ans = ans * x  +  *p++;
while( --i );

return( ans );
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl( x, coef, N )
double x;
double coef[];
int N;
{
double ans;
double *p;
int i;

p = coef;
ans = x + *p++;
i = N-1;

do
	ans = ans * x  + *p++;
while( --i );

return( ans );
}




















