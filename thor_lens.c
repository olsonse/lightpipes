/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h>
#include "pipes.h"

void lens();

void main(int argc, char *argv[]){
    void error_print();

    double F, x_shift, y_shift;

    /* Processing the command line argument  */
    if (argc<2 || argc >4){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&F);
    x_shift=y_shift=0.;
    if(argc == 3) sscanf(argv[2],"%le",&y_shift);
    if(argc == 4) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
    }
    read_field();

    if (field.double1 != 0.) {
        fprintf(stderr,"Lens can not be applied in spherical\
 coordinates,\nuse CONVERT first\n");
        exit(1);
    }

    
   
    lens(-F, x_shift, y_shift); 

    write_field();

}

void lens(F,xs,ys)

/* F>0 for negative lens !!!!! */

double F,xs,ys;
{ 
    int i,j,n2;
    long ik;
    double x,x2,y,dx,pi2, K;

    pi2=3.1415926*2.;
    K=pi2/field.lambda;
    n2=field.number/2;
    dx=field.size/field.number;

    ik=0;

    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx-xs;
        x2=x*x;
        for (j=1;j<=field.number; j++){ 
            double cab, sab, fi, cc;
            y=(j-n2-1)*dx-ys;
            fi=K*(x2+y*y)/(2.*F);
            cab=cos(fi);
            sab=sin(fi);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
            ik++;

        }

    }
}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s filters the beam through the quadratic phase corrector\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s  F[DX DY], where F is the focal length \n\
in [units you use],\n\
DX and DY are the transversal shifts of the lens optical axis \n\n",arr);

}






