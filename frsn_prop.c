/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

void fresnel(z)
double z;
{
int dims[2];
int i,j,fn2, fn22,io,jo,no2,dum,ii,ij,iiij,fresnl(), fftn();
long ik, ik1, ik2, ik3, ik4;
double *F_R, *F_I, *K_R, *K_I, RR, dx;
double cc, fc1, fs1, fc2, fs2, fc3, fs3, fc4, fs4, R1, R2, R3, R4;
double c4c1, c2c3, c4s1, s4c1, s2c3, c2s1, s4c3, s2c1, c4s3, s2s3, s2s1, c2s3, s4s1, c4c3, s4s3, c2c1, sh;

/* fprintf(stderr,"%e %d %e\n", z,n_old, dx); */
dx=field.size/(field.number-1.);
/*  Allocating a LOT OF MEMORY */

fn2=field.number*2;

F_R = (double *) calloc( fn2*fn2, sizeof(double) );
    if(F_R == NULL){fprintf(stderr,"Allocation error in Fresnel, \
use smaller grid\n");			
exit(1);}

F_I = (double *) calloc( fn2*fn2, sizeof(double) );
    if(F_I == NULL){fprintf(stderr,"Allocation error in Fresnel, \
use smaller grid\n");			
exit(1);}

K_R = (double *) calloc( fn2*fn2, sizeof(double) );
    if(K_R == NULL){fprintf(stderr,"Allocation error in Fresnel, \
use smaller grid\n");			
exit(1);}

K_I = (double *) calloc( fn2*fn2, sizeof(double) );
    if(K_I == NULL){fprintf(stderr,"Allocation error in Fresnel, \
use smaller grid\n");			
exit(1);}
/*****************************************************/

sh= +.5;
fn22=field.number+1;
no2=field.number/2;
RR=sqrt(1./(2*field.lambda*z))*dx*2;    

    ii=ij=1;
    ik=0;
    
    for (i=fn22-no2;i <= fn22+no2-1; i++){
      io=i-fn22;
      

 
      R1=RR*(io - .5 + sh);
      R3=RR*(io + .5 + sh);
      
 
       for (j=fn22-no2;j <= fn22+no2-1; j++){
	 iiij=ii*ij;
	  jo=j-fn22;
	  ik1=(i-1)*fn2+j-1;

/* Fresnel staff  */
	  R2=RR*(jo - .5 + sh);
	  R4=RR*(jo + .5 + sh);
	  dum=fresnl(R1,&fs1, &fc1);
	  dum=fresnl(R2,&fs2, &fc2);
	  dum=fresnl(R3,&fs3, &fc3);
	  dum=fresnl(R4,&fs4, &fc4);

	  c4c1=fc4*fc1;
	  c2s3=fc2*fs3;
	  c4s1=fc4*fs1;
	  s4c1=fs4*fc1;
	  s2c3=fs2*fc3;
	  c2s1=fc2*fs1;
	  s4c3=fs4*fc3;
	  s2c1=fs2*fc1;
	  c4s3=fc4*fs3;

	  s2s3=fs2*fs3;
	  s2s1=fs2*fs1;
	  c2c3=fc2*fc3;
	  s4s1=fs4*fs1;
	  c4c3=fc4*fc3;
	  c4c1=fc4*fc1;
	  s4s3=fs4*fs3;
	  c2c1=fc2*fc1;
	  

	  K_R[ik1]=0.5*(c4s3+s4c3-c4s1-s4c1-c2s3-s2c3+c2s1+s2c1)*iiij;
	  K_I[ik1]=0.5*(-c4c3+s4s3+c4c1-s4s1+c2c3-s2s3-c2c1+s2s1)*iiij;
/* Field staff */	
	  
/*
	 F_R[ik1]=(field.real[ik]-field.real[ik-1]+field.real[ik-field.number]-field.real[ik-field.number-1])*iiij*0.25;
	 F_I[ik1]=(field.imaginary[ik]-field.imaginary[ik-1]+field.imaginary[ik-field.number] - field.imaginary[ik-field.number-1])*0.25*iiij;
*/

F_R[ik1]=(field.real[ik])*iiij;
F_I[ik1]=(field.imaginary[ik])*iiij;
	
/*	  fprintf(stderr,"%ld\n",ik); */
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }
    dims[0]=dims[1]=fn2;
 
    fftn(2, dims, K_R, K_I, 1, (double) fn2);
  
    fftn(2, dims, F_R, F_I, 1, (double) fn2);
  
    ik=0;
    ii=ij=1;
    for(i=1; i<=fn2; i++){
     
      for(j=1; j<=fn2; j++){
	iiij=ii*ij;
	
	cc = K_R[ik]*F_R[ik]-K_I[ik]*F_I[ik];
	F_I[ik] = (K_R[ik]*F_I[ik]+F_R[ik]*K_I[ik])*iiij;
	F_R[ik]=cc*iiij;
	ik++;
	ij=-ij;
      } 
     ii=-ii;
    }
    
    free(K_R);
    free(K_I);


    fftn(2, dims, F_R, F_I, -1, 1.);
 
    ik=0;
    ii=ij=1;
    for(i=fn22-no2; i<=fn22+no2-1; i++){
      
      for(j=fn22-no2; j<=fn22+no2-1; j++){
	ik1=(i-1)*fn2+j-1;
	ik2=(i-2)*fn2+j-1;
	ik3=(i-2)*fn2+j-2;
	ik4=(i-1)*fn2+j-2;
	iiij=ii*ij;

	field.real[ik]=0.25*(F_R[ik1]-F_R[ik2]+F_R[ik3]-F_R[ik4])*iiij;
	field.imaginary[ik]=0.25*(F_I[ik1]-F_I[ik2]+F_I[ik3]-F_I[ik4])*iiij;  
	
/*	field.real[ik]=F_R[ik1]*iiij;
	field.imaginary[ik]=F_I[ik1]*iiij;
	*/	     
	ik++;
	ij=-ij;
      }
      ii=-ii;
    }
    free(F_R);
    free(F_I);
    return;
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




