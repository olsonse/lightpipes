/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"
#include <math.h>

void main(int argc, char *argv[]){
    void error_print();

    int i,j;
    double gain, i_sat, base, ss, intensity, dx,dx2;
    long ik1;



    /* Processing the command line argument  */

    if (argc != 4){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le", &gain);
    sscanf(argv[2],"%le",&base);
    sscanf(argv[3],"%le",&i_sat);  
   
    read_field();
    dx =field.size/(field.number);
    dx2 = dx*dx;
    
   
    ik1=0;
    for (i=1;i<=field.number ;i++){
        for (j=1;j<=field.number ;j++){
	  intensity=(field.real[ik1]*field.real[ik1]+\
field.imaginary[ik1]*field.imaginary[ik1]);


           ss=exp(base*(gain/(1.+(intensity/i_sat))));

                field.real[ik1] *= ss;
                field.imaginary[ik1] *= ss;
		   ik1++;
            }

      }
    write_field();


}


void error_print(char *arr) { 

    fprintf(stderr,"\n%s filters the field through a laser medium\n",arr );


    fprintf(stderr,"\nUSAGE: %s G L I, where G is the gain,\n\
L and I are the length and the saturation intensity of the medium\n",arr);

}












