/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void main(int argc, char *argv[]){

    double size_new, dx_new, dx_old, x_new, x_old, \
y_new, y_old, old_number, newfloat, inv_squares(), lower,\
upper, on21, nn21, x_shift, y_shift, angle, ss, cc, x0, y0, magnif, change, int16();
 long  n_old_max;

double z00, z01, z02 ,z03,\
       z10, z11 ,z12, z13,\
       z20 ,z21, z22, z23,\
       z30, z31, z32, z33;

int new_number, i, j, i_old, j_old;
void error_print();
   

#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);
#endif  

    /* Processing the command line argument  */

    if (argc< 2 || argc >7 ){
        error_print(argv[0]);
        exit(1);
    }
    read_field();
    /* reading  data from the command line */
    size_new=field.size;
    if((strstr(argv[1], "sam"))!= NULL ) {}
    else{
    if((sscanf(argv[1],"%le",&size_new))!=0){}
    else size_new=field.size;
    }


 
    old_number=field.number;
    n_old_max=old_number*(old_number-1)-1;
    new_number=field.number;
    if (argc >= 3){
  if((strstr(argv[2], "sam"))!= NULL ) {}
    else{ 
    if((sscanf(argv[2],"%d",&new_number))!=0){}
    else new_number=field.number;
    }
    }

    x_shift=y_shift=angle=0.;
    if(argc>=4) sscanf(argv[3],"%le",&y_shift);
    if(argc>=5) sscanf(argv[4],"%le",&x_shift);
    if(argc>=6) sscanf(argv[5],"%le",&angle);
    if(argc>=7) sscanf(argv[6],"%le",&magnif);
    else magnif=1.;
    if(magnif==0.){
        error_print();
        exit(1);
    }
    dx_new=size_new/(new_number-1.); 
    dx_old=field.size/(old_number-1.);
	change= 1.;
    angle *= Pi2/360.;
    field.number=new_number;
    field.size=size_new;

    if(fwrite ((char *)&(field.number), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.number\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.size), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.size\n");
        exit(1);
    }
    if(fwrite ((char *)&(field.lambda), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.lambda\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int1), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.int1\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int2), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.int2\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.int3), sizeof(int), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.int3\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double1), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.double1\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double2), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.double2\n");
        exit(1);
    }

    if(fwrite ((char *)&(field.double3), sizeof(double), 1, stdout) != 1){
        fprintf(stderr,"interpol: Error while writing FIELD.double3\n");
        exit(1);
    }

    on21= (int) old_number/2+1;
    nn21= (int) new_number/2+1;
    lower= (1-on21)*dx_old;
    upper= (old_number-on21)*dx_old;

    cc=cos(angle);
    ss=sin(angle);

    for (i=1; i<=new_number; i++){
        for (j=1; j<=new_number; j++){
	 
	  long n_tmp;
	  int i_small, i_local;

            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;

	    i_small=0;
	    i_local = 0;
/*          first row  */	   
	    n_tmp=(i_old-2)*old_number+j_old-2;

	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;} 
	    if (i_local != 1) z00=field.real[n_tmp];
	    else z00=0;
	    i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z10=field.real[n_tmp];
	    else z10=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z20=field.real[n_tmp];
	    else z20=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z30=field.real[n_tmp];
	    else z30=0;i_local=0;

/*          second row */
	    n_tmp=(i_old-2)*old_number+j_old-1;
	  if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z01=field.real[n_tmp];
	    else z01=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z11=field.real[n_tmp];
	    else z11=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z21=field.real[n_tmp];
	    else z21=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z31=field.real[n_tmp];
	    else z31=0;i_local=0;

/*          third row */
	    n_tmp=(i_old-2)*old_number+j_old;
	  if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z02=field.real[n_tmp];
	    else z02=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z12=field.real[n_tmp];
	    else z12=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z22=field.real[n_tmp];
	    else z22=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z32=field.real[n_tmp];
	    else z32=0;i_local=0;


/*          fourth     row */
	    n_tmp=(i_old-2)*old_number+j_old+1;
	    if( n_tmp < 0 || n_tmp > n_old_max){i_small=1; i_local=1;}
	    if (i_local != 1) z03=field.real[n_tmp];
	    else z03=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z13=field.real[n_tmp];
	    else z13=0;i_local=0;

            n_tmp += old_number;
	   if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z23=field.real[n_tmp];
	    else z23=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z33=field.real[n_tmp];
	    else z33=0;i_local=0;




	  if (i_small == 1){
	 
if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
newfloat=inv_squares(x_old,y_old,dx_old,z11,z21,z12,z22,x_new,y_new)/magnif;
            }
            else newfloat=0.;
           
	  }
	
else{ 

if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
newfloat=int16(\
x_old, y_old, x_new, y_new, dx_old,\
z00, z01, z02 ,z03,\
z10, z11 ,z12, z13,\
z20 ,z21, z22, z23,\
z30, z31, z32, z33)/magnif;}
else newfloat=0.;
	    
}
 if(fwrite (&newfloat, sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"interpol: error while writing FIELD.real  \n");
                exit (1);
        }
/* fprintf(stderr,"i= %d j= %d i_small=%d\n", i, j, i_small);	   */
	}
    }

/* phase is here */
   for (i=1; i<=new_number; i++){
        for (j=1; j<=new_number; j++){
	 
	  long n_tmp;
	  int i_small, i_local;

            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;

	    i_small=0;
	    i_local = 0;
/*          first row  */	   
	    n_tmp=(i_old-2)*old_number+j_old-2;

	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;} 
	    if (i_local != 1) z00=field.imaginary[n_tmp];
	    else z00=0;
	    i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z10=field.imaginary[n_tmp];
	    else z10=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z20=field.imaginary[n_tmp];
	    else z20=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z30=field.imaginary[n_tmp];
	    else z30=0;i_local=0;

/*          second row */
	    n_tmp=(i_old-2)*old_number+j_old-1;
	  if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z01=field.imaginary[n_tmp];
	    else z01=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z11=field.imaginary[n_tmp];
	    else z11=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z21=field.imaginary[n_tmp];
	    else z21=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z31=field.imaginary[n_tmp];
	    else z31=0;i_local=0;

/*          third row */
	    n_tmp=(i_old-2)*old_number+j_old;
	  if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z02=field.imaginary[n_tmp];
	    else z02=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z12=field.imaginary[n_tmp];
	    else z12=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z22=field.imaginary[n_tmp];
	    else z22=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z32=field.imaginary[n_tmp];
	    else z32=0;i_local=0;


/*          fourth     row */
	    n_tmp=(i_old-2)*old_number+j_old+1;
	    if( n_tmp < 0 || n_tmp > n_old_max){i_small=1; i_local=1;}
	    if (i_local != 1) z03=field.imaginary[n_tmp];
	    else z03=0;i_local=0;

            n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z13=field.imaginary[n_tmp];
	    else z13=0;i_local=0;

            n_tmp += old_number;
	   if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z23=field.imaginary[n_tmp];
	    else z23=0;i_local=0;

	    n_tmp += old_number;
	    if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
	    if (i_local != 1) z33=field.imaginary[n_tmp];
	    else z33=0;i_local=0;




	  if (i_small == 1){
	 
if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
newfloat=inv_squares(x_old,y_old,dx_old,z11,z21,z12,z22,x_new,y_new)/magnif;
            }
            else newfloat=0.;
           
	  }
	
else{ 
if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
newfloat=\
int16(\
x_old, y_old, x_new, y_new, dx_old,\
z00, z01, z02 ,z03,\
z10, z11 ,z12, z13,\
z20 ,z21, z22, z23,\
z30, z31, z32, z33)/magnif;}
 else newfloat=0.;	    
}
 if(fwrite (&newfloat, sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"interpol: error while writing FIELD.imaginary  \n");
                exit (1);
        }
/* fprintf(stderr,"i= %d j= %d i_small=%d\n", i, j, i_small);	  */
	}
    }



    }




void error_print(char *arr)
{
    fprintf(stderr,"\n%s  transfers the field into a grid with different size and dimension\n",arr);

    fprintf(stderr,"\nUSAGE: %s B [N, X_shift, Y_shift, A, M]  where \n\
B is the new size of the grid and  N is new dimension\n\
you can substitute the word <same> instead of B and N\n\
to preserve the existing values (if you only want to shift, rotate or scale)\n\
X_shift and Y_shift are the shifts of the field in a new grid\n\
A is the angle of rotation (after shifts) and\n\
M (M><0) is the magnification (the last applied)\n\n",arr);

}





double int16(\
xc, yc, xp, yp, dx,\
z00, z01, z02 ,z03,\
z10, z11 ,z12, z13,\
z20 ,z21, z22, z23,\
z30, z31, z32, z33)

double  xc, yc, xp, yp, dx,\
z00, z01, z02 ,z03,\
z10, z11 ,z12, z13,\
z20 ,z21, z22, z23,\
z30, z31, z32, z33;
{ double z1, z2, z3, z4, zz1, zz2, zz3, zz4, int4();

z1=z00; z2=z01; z3=z02, z4=z03;
zz1=int4(yc, dx, z1, z2, z3, z4, yp);

z1=z10; z2=z11; z3=z12, z4=z13;
zz2=int4(yc, dx, z1, z2, z3, z4, yp);

z1=z20; z2=z21; z3=z22, z4=z23;
zz3=int4(yc, dx, z1, z2, z3, z4, yp);

z1=z30; z2=z31; z3=z32, z4=z33;
zz4=int4(yc, dx, z1, z2, z3, z4, yp);

return int4(xc,dx,zz1,zz2,zz2,zz3,xp);

}

/* cubic four-point interpolation 

a+b*x1+c*x1^2+d*x1^3=y1
......................
a+b*x4+c*x4^2+d*x4^3=y4

where 
x1=x2-dx; 
x3=x2+dx; 
x4=x2+2*dx; 
the grid is uniform

*/

double int4(x2,dx,y1,y2,y3,y4,xz)
double dx,x2,y1,y2,y3,y4,xz;
{
double  t0,a,b,c,d;
 
/*    if(xz < x2 || xz >x2+dx ){
        fprintf(stderr,"out of range in the cubic interpolation routine \n");
        exit(1);
    }*/

/* maple staff, sorry */
t0=1./(dx*dx*dx); 

a = t0*(3.0*y2*x2*dx*dx+x2*x2*x2*y1+3.0*x2*x2*y1*dx+2.0*x2*y1*\
dx*dx-x2*x2*x2*y4+3.0*x2*x2*x2*y3+3.0*x2*x2*y3*dx-6.0*x2*y3*dx*dx+x2*y4*dx*dx+\
6.0*y2*dx*dx*dx-3.0*y2*x2*x2*x2-6.0*y2*x2*x2*dx)/6;

b = -t0*(-6.0*y3*dx*dx+y4*dx*dx+3.0*y2*dx*dx+9.0*x2*x2*y3-3.0*\
x2*x2*y4-9.0*y2*x2*x2-12.0*y2*x2*dx+6.0*y3*x2*dx+3.0*x2*x2*\
y1+6.0*y1*x2*dx+2.0*\
y1*dx*dx)/6;

d = -t0*(-3.0*y2-y4+y1+3.0*y3)/6;

c = t0*(-3.0*y2*x2+y3*dx+y1*dx-2.0*y2*dx-x2*y4+3.0*x2*y3+x2*y1)/2;


 
/*return a+b*xz+c*xz*xz+d*xz*xz*xz;*/
 return a+xz*(b+xz*(c+xz*d)); 

}

double inv_squares(x,y,dx,z,zx,zy,zxy,x1,y1)
double x,y,dx,z,zx,zy,zxy,x1,y1;
{
    double tol;
    double s1,s2,s3,s4,xlow,xhigh,ylow,yhigh, sum;
    tol=1e-6*dx;
    if(x1< x-tol || x1>x+dx+tol || y1<y-tol || y1>y+dx+tol){
        fprintf(stderr,"out of range in inv_squares %g %g %g %g %g %g %g\n", x,x1,y,y1,dx, y1-y, x1-x);
        exit(1);
    }

    xlow=x1-x;
    xhigh=x+dx-x1;
    ylow=y1-y;
    yhigh=y+dx-y1;

    if(xlow< -tol || xhigh< -tol || ylow < -tol || yhigh < -tol ){ 
        fprintf(stderr," inv_squares: out of range, %g %g %g %g\n", xlow,xhigh, ylow, yhigh);
        exit(1);
    }


    if (fabs(xlow) < tol) return z+ylow*(zy-z)/dx;
    if (fabs(ylow) < tol) return z+xlow*(zx-z)/dx;
    if (fabs(xhigh) < tol) return zx+ylow*(zxy-zx)/dx;
    if (fabs(yhigh) < tol) return zy+xlow*(zxy-zx)/dx;

    s1=1./(xlow*ylow);
    s2=1./(xhigh*ylow);
    s3=1./(xlow*yhigh);
    s4=1./(xhigh*yhigh);


    sum=s1+s2+s3+s4;
    s1=s1/sum;
    s2=s2/sum;
    s3=s3/sum;
    s4=s4/sum;

    return z*s1+zx*s2+zy*s3+zxy*s4;

}





















