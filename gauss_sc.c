/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h>
#include "pipes.h"

void gauss_scr();

void main(int argc, char *argv[]){
    void error_print();
    double R, x_shift, y_shift, AA;

    /* Processing the command line argument  */
    if (argc<2 || argc >5){
        error_print();
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&R);
    x_shift=y_shift=0.;
    AA=1.;
    if(argc == 3) sscanf(argv[2],"%le",&y_shift);
    if(argc == 4) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
    }
    if(argc == 5) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
	sscanf(argv[4],"%le",&AA);
    }
    read_field();

    gauss_scr(R, x_shift, y_shift, AA); 

    write_field();

}

void gauss_scr(R,xs,ys, AA)
double R,xs,ys, AA;
/* 1-AA is the transmission in the minimum point 
   if the screen is modeled as a mirror then AA
   is the maximum reflection and gauss_scr returns the
   transmitted field
*/
{ 
    int i,j,n2;
    long ik;
    double x,x2,y,y2,dx,cc,R2;

    n2=field.number/2;
    dx=field.size/field.number;
    R2=R*R;
      
    ik=0;

    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx-xs;
        x2=x*x;
        for (j=1;j<=field.number; j++){
            y=(j-n2-1)*dx-ys;
            y2=y*y;
            cc=sqrt(fabs(1.- AA*exp(-(x2+y2)/R2)));
            field.imaginary[ik] *= cc;
            field.real[ik] *= cc;
            ik++;
        }

    }
}




void error_print()
{
    fprintf(stderr,"\ngauss_screen filters the beam through a gaussian screen\n");
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"gauss_screen R [DX DY T], where R is the 1/e (intensity) \n\
radius of the screen in [units you use],\n\
DX and DY are the shifts of the screen  center\n\
(1-T) gives the minimum transmission, default T=1\n\
i.e the default  minimum transmission == 0\n\n");

}





