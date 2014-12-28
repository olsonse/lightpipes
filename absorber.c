/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"
#include "math.h"

void main(int argc, char *argv[]){
    void error_print();

    int i,j;
    double R, rr;
    long ik1;


    /* Processing the command line argument  */

    if (argc<2 || argc >2){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&R);

    rr=sqrt(R);
    read_field();
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            ik1=(i-1)*field.number+j-1; 

            field.real[ik1] *= rr;
            field.imaginary[ik1] *= rr; 


        }
    }


    write_field();


}


void error_print(char *arr)
{
    fprintf(stderr,"\n%s filters the field through \
intensity attenuator R\n", arr);

    fprintf(stderr,"\nUSAGE: %s K, where K \
is the coefficient \nof attenuation or amplification (must be positive)\n\n", arr);


}


