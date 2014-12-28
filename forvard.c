/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
#include <math.h>
#include "pipes.h"
#include "fftn.h"
#undef REAL
#define REAL double

void fft3();
int fftn ();

void main(int argc, char *argv[]){
    void error_print();
    void forvard();

    double z;

    /* Processing the command line argument  */
    if (argc!=2){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&z);
    read_field();
    if (field.double1 !=0.) {
        fprintf(stderr,"Forvard can not be applied in spherical\
 coordinates,\nuse CONVERT first\n");
        exit(1);
    }

    forvard(z); 




    write_field();


}

void error_print(char *arr)
{
    fprintf(stderr,"\n%s propagates the field to \
distance Z [units you use]\n\
using FFT algorithm",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s Z, where  Z is the distance to propagate\n\n",arr);

}


#include "fft_prop.c"
