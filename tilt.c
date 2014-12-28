/*--------------------------------------------------------------*/
/*	(C) Gleb Vdovin 1993-1995				*/
/*	This file is a part of LightPipes beta 			*/
/*	(October 1995) distribution				*/
/*	Send bug reports to gleb@okotech.com     		*/
/*								*/
/*	This file may be distributed ONLY together with		*/
/*	the complete source of LightPipes  			*/	
/*--------------------------------------------------------------*/

#include <math.h>
#include "pipes.h"

void tilt();

void main(int argc, char *argv[]){
    void error_print();

    double tx,ty;

    /* Processing the command line argument  */
    if (argc!=3){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
   
    sscanf(argv[1],"%le",&ty);
    sscanf(argv[2],"%le",&tx);
    read_field();

    tilt(tx,ty);
/*    fprintf(stderr,"%g %g ", tx,ty);*/

    write_field();

}

void tilt(tx,ty)
double tx,ty;
{ 
    int i,j,n2;
    long ik;
    double x,y,dx,fi,K,cab,sab,cc;

    n2=field.number/2;
    K=2.*3.141592654/field.lambda;
    dx=field.size/field.number;

    ik=0;
    /*     fprintf(stderr,"%le %le\n", tx,ty);*/ 
    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx;
        for (j=1;j<=field.number; j++){
            y=(j-n2-1)*dx;
            fi= -(tx*x+ty*y)*K;
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
    fprintf(stderr,"\n%s introduces tilt into the field distribution\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s X Y, where  X and Y  are the X and Y components in \
radians\n\n",arr);

}




