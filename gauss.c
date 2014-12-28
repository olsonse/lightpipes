/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h>
#include "pipes.h"

void gauss();

void main(int argc, char *argv[]){
    void error_print();
    double R, x_shift, y_shift, A;

    /* Processing the command line argument  */
    if (argc<2 || argc >5){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&R);
    x_shift=y_shift=0.;
    A=1.;
    if(argc == 3) sscanf(argv[2],"%le",&y_shift);
    if(argc == 4) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
    }
    if(argc == 5) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
	sscanf(argv[4],"%le",&A);
    }
    read_field();

    gauss(R, x_shift, y_shift, A); 

    write_field();

}

void gauss(R,xs,ys,AA)
double R,xs,ys, AA;
{ 
    int i,j,n2;
    long ik;
    double x,x2,y,y2,dx,cc,R2;

    n2=field.number/2;
    dx=field.size/field.number;
    R2=R*R*2.;
    AA=sqrt(fabs(AA));
    ik=0;

    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx-xs;
        x2=x*x;
        for (j=1;j<=field.number; j++){
            y=(j-n2-1)*dx-ys;
            y2=y*y;
            cc=AA*exp(-(x2+y2)/R2);
            field.imaginary[ik] *= cc;
            field.real[ik] *= cc;
            ik++;
        }

    }
}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s filters the beam through a gaussian diaphragm\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s R [DX DY T], where R is the 1/e (intensity) \n\
radius of the diaphragm in [units you use],\n\
DX and DY are the shifts of the diaphragm  center\n\
T is the intensity transmission in the  maximum (default T=1)\n\n",arr);

}




