/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"

void main(int argc, char *argv[]){
    void error_print();

    int i,j,i2;
    double R, dx,x,y,rr, x_shift, y_shift;
    long ik1;


    /* Processing the command line argument  */

    if (argc<2 || argc >4){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&R);
    x_shift=y_shift=0.;
    if(argc == 3) sscanf(argv[2],"%le",&y_shift);
    if(argc == 4) {
        sscanf(argv[2],"%le",&y_shift);
        sscanf(argv[3],"%le",&x_shift);
    }
    rr=R*R;
    /*     fprintf(stderr,"%le\n",rr); */
    read_field();

    dx =field.size/(field.number);
    i2=field.number/2;

    /* Cuttitng the aperture      */


    for (i=1;i<=field.number; i++){
        x=(i-i2-1)*dx-x_shift;
        for (j=1;j<=field.number; j++){

            y=(j-i2-1)*dx-y_shift;
            ik1=(i-1)*field.number+j-1; 
            /*        printf("%lf %lf    ", x*x+y*y, rr); */
            if(x*x+y*y <= rr) {
                field.real[ik1]=0.;
                field.imaginary[ik1]=0.; 
            }
            /*           printf("%ld\n",  ik1); */


        }
    }


    write_field();


}


void error_print(char *arr)
{
    fprintf(stderr,"\n%s filters the field through \
nontransparent circular screen of radius R\n", arr);

    fprintf(stderr,"\nUSAGE: %s R [DX DY], where R \
is the radius \nof the screen in [units you use],\n\
Dx and DY are zhe shifts of the screen center\n\n",arr);


}


