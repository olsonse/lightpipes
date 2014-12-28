/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h>
#include "pipes.h"



void main(int argc, char *argv[]){
    void error_print();
 
    int i,j;
    long ik;
    double z, z1 , f, f1, ampl_scale;
    void fresnel();
    /* Processing the command line argument  */
    if (argc!=3){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&f);
    sscanf(argv[2],"%le",&z);
    
    if(f==0)fprintf(stderr,"lens_forvard: Lens with \
ZERO focal length does not exist.\n");

    read_field();


    f1=0.;
    if (field.double1 !=0. ) f1=1./field.double1;
    else f1=10000000.* field.size*field.size/field.lambda;
    if( (f+f1) != 0.) f=(f*f1)/(f+f1);
    else f=10000000.* field.size*field.size/field.lambda;
    z1=-z*f/(z-f);
    
if(z1 < 0. ) {fprintf(stderr,"Sorry, lens_fresn can not propagate behind\n\
the focal point, use lens_forvard instead, exiting.\n\n");
exit (1);
    }


    fresnel(z1); 


    ampl_scale=(f-z)/f;
    field.size *= ampl_scale;
    field.double1= -1./(z-f);
    

    ik=0;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]/ampl_scale;
            field.imaginary[ik] =field.imaginary[ik]/ampl_scale; 
            ik++;

        }

    }
    

  
 
    write_field();


}



void error_print(char *arr)
{
    fprintf(stderr,"\n%s propagates the field to \
distance Z [units you use]\nin VARIABLE COORDINATE  system",arr);
    fprintf(stderr,"\n\nUSAGE:  ");
    fprintf(stderr,"%s F Z, where  F is \
the focal lenght of the input lens,\n\
Z is the distance to propagate\n\n",arr);

}

#include "frsn_prop.c"


