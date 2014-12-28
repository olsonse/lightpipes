/*--------------------------------------------------------------*/
/*	(C) Gleb Vdovin 1993-1999				*/
/*	Send bug reports to gleb@okotech.com     		*/
/*								*/
/*	This file may be distributed ONLY together with		*/
/*	the complete source of LightPipes  			*/	
/*--------------------------------------------------------------*/

#include "pipes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void main(int argc, char *argv[]){

    double size_new, dx_new, dx_old, x_new, x_old, \
y_new, y_old, old_number, newfloat, inv_squares(), lower,\
upper, on21, nn21, x_shift, y_shift, angle, ss, cc, x0, y0, magnif, change;
    int new_number,i,j,i_old, j_old;
    void error_print();
    long n_old, n_oldx, n_oldy, n_oldxy;
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
    /*    n_old_max=old_number*(old_number-1)-1;*/
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
            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;
            n_old=(i_old-1)*old_number+j_old-1;
            n_oldx=n_old+old_number;
            n_oldy=n_old+1;
            n_oldxy=n_oldx+1;
	   /*  fprintf(stderr,"%d %d\n", i_old, j_old);*/
            if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
                /*   if(x_new <x_old || x_new > x_old+dx_old) fprintf(stderr,"errorx\n");
                	    if(y_new <y_old || y_new > y_old+dx_old) fprintf(stderr,"errory %i %i \n",i,j);*/
                /*    fprintf(stderr,"%i %g %g %g \n", j, y_old, y_new, y_old+dx_old);*/

newfloat=inv_squares(x_old,y_old,dx_old,field.real[n_old],field.\
real[n_oldx],field.real[n_oldy],field.real[n_oldxy],x_new,y_new)/magnif;
newfloat *= change;
            }
            else newfloat=0.;
            if(fwrite (&newfloat, sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"interpol: error while writing FIELD.real  \n");
                exit (1);
            }

        }
    }


    for (i=1; i<=new_number; i++){
        for (j=1; j<=new_number; j++){
            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int) floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;
            n_old=(i_old-1)*old_number+j_old-1;
            n_oldx=n_old+old_number;
            n_oldy=n_old+1;
            n_oldxy=n_oldx+1;

            if(x_new > lower && x_new < upper && y_new >lower && y_new < upper){
                /*    if(x_new <x_old || x_new > x_old+dx_old) fprintf(stderr,"errorx1\n");
                	    if(y_new <y_old || y_new > y_old+dx_old) fprintf(stderr,"errory1 %i %i \n",i, j); */

                newfloat=inv_squares(x_old,y_old,dx_old,field.imaginary[n_old],field.\
imaginary[n_oldx],field.imaginary[n_oldy],field.imaginary[n_oldxy],\
x_new,y_new)/magnif;
newfloat *= change;
            }
            else newfloat=0.;
            if(fwrite (&newfloat, sizeof(double), 1, stdout)!=1){
                fprintf(stderr,"interpol: error while writing FIELD.real  \n");
                exit (1);
            }

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



/*
*         Inverse square interpolation : 
*         given square (x,y) (x+dx,y) (x,y+dx) (x+dx,y+dx) 
*         with values   z     zx       zy        zxy 
*         the program returns value for Z for arbitrary 
*         x1 and y1 inside the rectangle.
*/


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













