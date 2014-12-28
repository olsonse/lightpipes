/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
#include <math.h>
#include "pipes.h"
typedef struct {
    double r,i;
}
fcomplex;
fcomplex Cadd();
fcomplex Csub();
fcomplex Cmul();
fcomplex Complex();
fcomplex Conjg();
fcomplex Cdiv();
double Cabs();
fcomplex Csqrt();
fcomplex RCmul();

void all_mem(), free_mem(), error_print(), elim();

fcomplex * a, * b, * c, * alpha, * beta, * u, * p, *u1, *u2;
double *refr, *absorb;

void main(int argc, char *argv[]){
    void error_print();
    void step();
    long ik, ikij,  ikij1,  ikij_1, iki1j, iki_1j;
    double z, delta, delta2, Pi4lz,A,AA,band_pow;
    fcomplex uij, uij1, uij_1, ui1j, ui_1j, medium;
    int i,j,jj, ii, nstep, istep, i_left, i_right;
    FILE *fr;

   
       /* Processing the command line argument  */
    if (argc< 2 || argc >5 ){
        error_print();
        exit(1);
    }

 read_field();
    all_mem();

    /* reading the data from a command line */
    sscanf(argv[1],"%le",&z);
    nstep = 1;
    if (argc > 2 ) sscanf(argv[2],"%d", &nstep);



if (argc >3 ) {
if((strstr(argv[3], "void")) != NULL ) {
  fprintf(stderr,"Steps: void refractive coefficient file, skipping \n");
  goto l1;}
 
if((fr=fopen(argv[3],"r"))==NULL){
  fprintf(stderr,"Steps: error opening file %s, exiting \n", argv[3]);
  exit(1);} 

    ik=0;  
    for (i=1;i<= field.number; i ++){
        for (j=1;j <= field.number;j ++){
	double fi;
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"Steps: reading the refractive indices:  end of input file reached, exiting\n");
	      exit(1);}
	      
	      refr[ik]=fi;
	      ik ++;
	    	  
	}      
    }
fclose(fr);
}

l1: 
if (argc >4 ) {
  if((fr=fopen(argv[4],"r"))==NULL){
  fprintf(stderr,"Steps: error opening file %s, exiting \n", argv[4]);
  exit(1);} 

    ik=0;  
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
	double fi;
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"Steps: absorbtion  end of input file reached, exiting\n");   
	      exit(1);
	      }
	      absorb[ik]=fi;
	      ik ++;
	    	  
	}
    }

fclose(fr);
}

/* the arrays for refraction and absorption are finally formed */




    if (field.double1 !=0.) {
        fprintf(stderr,"Steps can not be applied in spherical\
 coordinates,\nuse CONVERT first\n");
        exit(1);
    }

    
    z=z/2.;
    Pi4lz = 4.*Pi/field.lambda/z;
    delta=field.size/((double)(field.number-1.));
    delta2 = delta*delta;
    A=0.;
/* absorbtion at the borders is described here */
    AA=10./z/nstep; /* total absorbtion */
    band_pow=2.;   /* profile of the absorbtion border, 2=quadratic*/
/* width of the absorbtion border */
    i_left=field.number/2+1-0.4*field.number;
    i_right=field.number/2+1+0.4*field.number;
/* end absorbtion */
  

    for (i=1; i<=field.number; i++){
      u2[i].r=u2[i].i = 0.;
/*      a[i].r=b[i].r= -1./delta2; */
      a[i].i=b[i].i=0.;
    }


    
/*  Main  loop, steps here */
    for(istep = 1; istep <= nstep ; istep ++){

/*  Elimination in the direction i, halfstep  */
    for(jj=2; jj<= field.number-2; jj += 2){

        j=jj;
        for (i=2; i<=field.number-1; i++){
	  
            ikij=(j-1)*field.number+i-1;
            ikij1=(j)*field.number+i-1;
            ikij_1=(j-2)*field.number+i-1;

    
            uij=Complex(field.real[ikij], field.imaginary[ikij]);
            uij1=Complex(field.real[ikij1], field.imaginary[ikij1]);
            uij_1=Complex(field.real[ikij_1], field.imaginary[ikij_1]);

            p[i] = RCmul(-2., uij);
            p[i] = Cadd(p[i], uij1);
            p[i] = Cadd(p[i], uij_1);
            p[i] = RCmul(-1./delta2/refr[ikij], p[i]);
            p[i] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[i]);

        }

	for (i=1; i<=field.number; i++){
	  ik=(j-1)*field.number + i -1;

/* 	  fprintf (stderr,"%d %d %e %e\n", j , i, refr[ik], absorb[ik]);  */

/*	  medium.r = -4*Pi*Pi*(refr[ik]*refr[ik]-1.)/field.lambda/field.lambda; */       medium.r =0.;
	  medium.i= 2.*Pi*absorb[ik]/field.lambda;
	  c[i].r=-2./delta2/refr[ik];
	  c[i].i= Pi4lz;
	  c[i].r += medium.r;
	  c[i].i += medium.i;  
	  a[i].r=b[i].r= -1./delta2/refr[ik];
/* absorbtion borders are formed here */
      if( i <= i_left){
int iii=i_left-i;
c[i].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(i_left)),band_pow);} 
      if( i >= i_right){ 
	int iii =i-i_right;
	int im=field.number-i_right;
c[i].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(im)),band_pow);}
/* end absorbtion */

}
        elim();

        for (i=1; i<=field.number; i++){
            ikij_1=(j-2)*field.number+i-1;
            field.real[ikij_1]=u2[i].r;
            field.imaginary[ikij_1]=u2[i].i;
            u1[i].r=u[i].r;
            u1[i].i=u[i].i;
        }

        j=jj+1;

        for (i=2; i<=field.number-1; i++){

            ikij=(j-1)*field.number+i-1;
            ikij1=(j)*field.number+i-1;
            ikij_1=(j-2)*field.number+i-1;


            uij=Complex(field.real[ikij], field.imaginary[ikij]);
            uij1=Complex(field.real[ikij1], field.imaginary[ikij1]);
            uij_1=Complex(field.real[ikij_1], field.imaginary[ikij_1]);

            p[i] = RCmul(-2., uij);
            p[i] = Cadd(p[i], uij1);
            p[i] = Cadd(p[i], uij_1);
            p[i] = RCmul(-1./delta2/refr[ikij], p[i]);
            p[i] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[i]);

        }
	for (i=1; i<=field.number; i++){
	  ik=(j-1)*field.number + i -1;

/* 	  medium.r = -4*Pi*Pi*(refr[ik]*refr[ik]-1.)/field.lambda/field.lambda;       */    medium.r =0.;
	  medium.i= 2.*Pi*absorb[ik]/field.lambda;
	  c[i].r=-2./delta2/refr[ik];
	  c[i].i= Pi4lz;
	  c[i].r += medium.r;
	  c[i].i += medium.i;   
a[i].r=b[i].r= -1./delta2/refr[ik];
/* absorbtion borders are formed here */
      if( i <= i_left){
int iii=i_left-i;
c[i].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(i_left)),band_pow);} 
      if( i >= i_right){ 
	int iii =i-i_right;
	int im=field.number-i_right;
c[i].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(im)),band_pow);}
/* end absorbtion */

}
        elim();

        for (i=1; i<=field.number; i++){
            ikij=(j-2)*field.number+i-1;
            field.real[ikij]=u1[i].r;
            field.imaginary[ikij]=u1[i].i;

            u2[i].r=u[i].r;
            u2[i].i=u[i].i;
        }
    }

    for (i=1; i<=field.number; i++){
        ikij=(field.number-1)*field.number+i-1;
        field.real[ikij]=u2[i].r;
        field.imaginary[ikij]=u2[i].i;

    } 

/* Elimination in the j direction is here, halfstep */

  for (i=1; i<=field.number; i++){
    u2[i].r=u2[i].i = 0.;}

    for(ii=2; ii<= field.number-2; ii += 2){

        i=ii;
        for (j=2; j<=field.number-1; j++){

            ikij=(j-1)*field.number+i-1;
            iki1j=(j-1)*field.number+i;
            iki_1j=(j-1)*field.number+i-2;

    
            uij=Complex(field.real[ikij], field.imaginary[ikij]);
            ui1j=Complex(field.real[iki1j], field.imaginary[iki1j]);
            ui_1j=Complex(field.real[iki_1j], field.imaginary[iki_1j]);

            p[j] = RCmul(-2., uij);
            p[j] = Cadd(p[j], ui1j);
            p[j] = Cadd(p[j], ui_1j);
            p[j] = RCmul(-1./delta2/refr[ikij], p[j]);
            p[j] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[j]);

        }
	for (j=1; j<=field.number; j++){
	  ik=(j-1)*field.number + i -1;

/*	  medium.r = -4*Pi*Pi*(refr[ik]*refr[ik]-1.)/field.lambda/field.lambda;    */    medium.r =0.;
	  medium.i= 2.*Pi*absorb[ik]/field.lambda;
	  c[j].r=-2./delta2/refr[ik];
	  c[j].i= Pi4lz;
	  c[j].r += medium.r;
	  c[j].i += medium.i;  
a[j].r=b[j].r= -1./delta2/refr[ik];
/* absorbtion borders are formed here */
      if( j <= i_left){
int iii=i_left-j;
c[j].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(i_left)),band_pow);} 
      if( j >= i_right){ 
	int iii =j-i_right;
	int im=field.number-i_right;
c[j].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(im)),band_pow);}
/* end absorbtion */

}
        elim();

        for (j=1; j<=field.number; j++){
            iki_1j=(j-1)*field.number+i-2;
            field.real[iki_1j]=u2[j].r;
            field.imaginary[iki_1j]=u2[j].i;
            u1[j].r=u[j].r;
            u1[j].i=u[j].i;
        }

        i=ii+1;
   for (j=2; j<=field.number-1; j++){

            ikij=(j-1)*field.number+i-1;
            iki1j=(j-1)*field.number+i;
            iki_1j=(j-1)*field.number+i-2;

    
            uij=Complex(field.real[ikij], field.imaginary[ikij]);
            ui1j=Complex(field.real[iki1j], field.imaginary[iki1j]);
            ui_1j=Complex(field.real[iki_1j], field.imaginary[iki_1j]);

            p[j] = RCmul(-2., uij);
            p[j] = Cadd(p[j], ui1j);
            p[j] = Cadd(p[j], ui_1j);
            p[j] = RCmul(-1./delta2/refr[ikij], p[j]);
            p[j] = Cadd(Cmul(Complex(0., Pi4lz), uij),p[j]);

        }
	for (j=1; j<=field.number; j++){
	  ik=(j-1)*field.number + i -1;

/*	  medium.r = -4*Pi*Pi*(refr[ik]*refr[ik]-1.)/field.lambda/field.lambda;     */    medium.r =0.;
	  medium.i= 2.*Pi*absorb[ik]/field.lambda;
	  c[j].r=-2./delta2/refr[ik];
	  c[j].i= Pi4lz;
	  c[j].r += medium.r;
	  c[j].i += medium.i;  
a[j].r=b[j].r= -1./delta2/refr[ik];
/* absorbtion borders are formed here */
      if( j <= i_left){
int iii=i_left-j;
c[j].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(i_left)),band_pow);} 
      if( j >= i_right){ 
	int iii =j-i_right;
	int im=field.number-i_right;
c[j].i += (AA*2.*Pi/field.lambda)*pow((double) iii/ ((double)(im)),band_pow);}
/* end absorbtion */

}
        elim();
       

        for (j=1; j<=field.number; j++){
            ikij=(j-1)*field.number+i-2;
            field.real[ikij]=u1[j].r;
            field.imaginary[ikij]=u1[j].i;

            u2[j].r=u[j].r;
            u2[j].i=u[j].i;
        }
    }

    for (j=1; j<=field.number; j++){
        ikij=(j-2)*field.number+i-1;
        field.real[ikij]=u2[j].r;
        field.imaginary[ikij]=u2[j].i;

    }
    
/* end j */

    }
 write_field();

    free_mem();
    
    
}


void elim()
{ 
    int i;
    fcomplex cc;
    /* initial condition, everything is going to be zero at the edge */
    alpha[2].r=0.;
    alpha[2].i=beta[2].r=beta[2].i=0.;

    alpha[field.number].r=0.;
    alpha[field.number].i=beta[field.number].r=beta[field.number].i=0.;

    /* forward elimination */
    for(i=2;i<=field.number-2; i++){
        cc=Csub(c[i],Cmul(a[i],alpha[i]));
        alpha[i+1]=Cdiv(b[i],cc);
        beta[i+1]=Cdiv(Cadd(p[i],Cmul(a[i],beta[i])),cc);
    }
    i=field.number;
    cc=Csub(c[i],Cmul(a[i],alpha[i]));
    beta[field.number+1]=Cdiv(Cadd(p[i],Cmul(a[i],beta[i])),cc);
    /* edge amplitude =0 */
    u[field.number]=beta[field.number+1];
    /* backward elimination        */
    for(i=field.number-1; i>=1; i--){
        u[i]=Cadd(Cmul(alpha[i+1],u[i+1]),beta[i+1]);

    }
}

void all_mem()
{ int i,j; long ik;
/*    Allocatiom of the  memory for arrays: */

if (( refr=(double *) calloc( (field.number+2)*(field.number+2), sizeof(double) ))==NULL){
        fprintf(stderr,"Steps: allocation error, array of refractive indices");
        exit(1);
    }

if (( absorb=(double *) calloc( (field.number+2)*(field.number+2), sizeof(double) ))==NULL){
        fprintf(stderr,"Steps allocation error, array of absorbtion");
        exit(1);
    }

ik=0;
    for(i=1; i<= field.number+2; i++){
      for(j=1; j<= field.number+2; j++){
	refr[ik]=1.;
	absorb[ik]=0.;
	ik ++;
      }
}

/* end arrays */
    /* allocation of memory for the elimination */

    if (( a=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    if (( p=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    if (( b=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }
    if (( c=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    if (( u=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    if (( u1=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    if (( u2=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }
    if (( alpha=(fcomplex *) calloc( field.number+3, sizeof(fcomplex)))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }
    if (( beta=(fcomplex *) calloc( field.number+3, sizeof(fcomplex) ))==NULL){
        fprintf(stderr,"Allocation error");
        exit(1);
    }

    /* Memory done */ 

}

void free_mem()

{
    /* freeing  the memory        */

    free(a);
    free(b);
    free(c);
    free(u);
    free(u1);
    free(u2);
    free(alpha);
    free(beta);
    free(p);
    free(absorb);
    free(refr);

}


/* traditional - K&R */

fcomplex Cadd(a,b)
fcomplex a,b;
{
    fcomplex c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplex Csub(a,b)
fcomplex a,b;
{
    fcomplex c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplex Cmul(a,b)
fcomplex a,b;
{
    fcomplex c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplex Complex(re,im)
double im,re;
{
    fcomplex c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplex Conjg(z)
fcomplex z;
{
    fcomplex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplex Cdiv(a,b)
fcomplex a,b;
{
    fcomplex c;
    double den;

    den=b.r*b.r+b.i*b.i;

    if( den != 0.){
        c.r= (a.r*b.r +a.i*b.i)/den;
        c.i= (a.i*b.r - a.r*b.i)/den;
    }
    else{
        fprintf(stderr,"Complex division by zero, exiting\n");
        exit(1);
    }

    return c;
}

double Cabs(z)
fcomplex z;
{
    double x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    }
    else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplex Csqrt(z)
fcomplex z;
{
    fcomplex c;
    double x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    }
    else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        }
        else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        }
        else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplex RCmul(x,a)
fcomplex a;
double x;
{
    fcomplex c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}

void error_print() { 

fprintf(stderr,"\nsteps propagates the field to\n\
distance N*Z [units you use]\n\
using finite difference method\n");

fprintf(stderr,"\nUSAGE: \n"); 

fprintf(stderr,"steps Z [N R A] ,\n\
where Z is the step size, N is the number of steps to perform\n\
R is the name of the file, which  contains the distribution\n\
of refracive index in gnuplot format\n\
if R equals to < void > then this option is skipped\n\
A is the name of a file which  contains the two-dimensional\n\
distribution of the absorbtion coefficient\n\n");

fprintf(stderr,"Examples:\n\
steps 0.1 20 - performs 20 steps of 0.1 m in a free space,\n\
total distance of propagation is  2 m \n\n\
steps 0.1 20 void abs_d - performs 20 steps of 0.1 m in a medium\n\
absorbtion coefficient of which is given in a file abs_d,\n\
refraction coefficient equals to unity\n\n");


fprintf(stderr,"steps 0.1 20 refr_d abs_d - performs 20 steps of 0.1 m in a medium\n\
absorbtion coefficient of which is given in a file abs_d,\n\
refraction coefficient is given in a file refr_d\n\n\
steps 0.1 20 refr_d - performs 20 steps of 0.1 m in a medium\n\
refraction coefficient of which is given in a file refr_d\n\
without absorbtion\n\n");

}








