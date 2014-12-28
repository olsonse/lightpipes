/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"
#include <math.h>
#include <string.h>

void main(int argc, char *argv[]){
    void error_print();

    int i,j;
    double sum, dx, dx2, asum;
    long ik1;



    /* Processing the command line argument  */

    if (argc!= 2){
        error_print(argv[0]);
        exit(1);
    }


    read_field();
    sum=0;
    dx =field.size/(field.number);
    dx2 = dx*dx;



    /* Calculating the power */

    for (i=1;i<=field.number ;i++){
        for (j=1;j<=field.number ;j++){

            ik1=(i-1)*field.number+j-1;

            sum += (field.real[ik1]*field.real[ik1]+ \
		    field.imaginary[ik1]*field.imaginary[ik1]) /* *dx2*/ ;

        }
}
    sum *= dx2;

    if (sum == 0) {
        fprintf(stderr,"normal: Zero beam power, program terminated\n");
        exit(1);
    }
    asum=sqrt(1./sum);
    /*Normalizing the power */
    for (i=1;i<=field.number ;i++)
        for (j=1;j<=field.number ;j++){

            ik1=(i-1)*field.number+j-1;

            field.real[ik1] *= asum;
            field.imaginary[ik1] *=asum;

        }
    /*    sscanf(argv[1],"%c", &in);*/
    if (strstr(argv[1], "y")!= NULL)fprintf(stderr,"normal: Normalization coefficient= %e\n",sum);

    write_field();


}


void error_print( char *arr) { 

    fprintf(stderr,"\n%s normalizes the beam power to unity and\
\noptionally prints the normalisation coefficient\n \n",arr);


    fprintf(stderr,"\nUSAGE: %s [y,n] , if you pass y \n\
then the normalization coefficient is printed to stderr \n\
if you pass n - normalization is silent\n\n",arr);




}
