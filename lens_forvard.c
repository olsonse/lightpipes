/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    


#include <math.h>
#include "pipes.h"

void fft3();
void fftn ();

void main(int argc, char *argv[]){
    void error_print();
    void forvard();
    int i,j;
    long ik;
    double z, z1 , f, f1, ampl_scale, *fi1, *fr1;

    /* Processing the command line argument  */
    if (argc!=3){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&f);
    sscanf(argv[2],"%le",&z);
    
  
  if(f==0){ fprintf(stderr,"lens_forvard: Lens with \
ZERO focal length DOES NOT EXIST.\n");
	      exit (1);
	    }
    read_field();


    f1=0.;
    if (field.double1 !=0. ) f1=1./field.double1;
    else f1=10000000.* field.size*field.size/field.lambda;
    if( (f+f1) != 0.) f=(f*f1)/(f+f1);
    else f=10000000.* field.size*field.size/field.lambda;
    z1=-z*f/(z-f);


    forvard(z1); 


    ampl_scale=(f-z)/f;
    field.size *= ampl_scale;
    field.double1= -1./(z-f);

    if (z1 >= 0.){
    ik=0;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]/ampl_scale;
            field.imaginary[ik] =field.imaginary[ik]/ampl_scale; 
            ik++;

        }

    }
    }
    else{
      fr1 = (double *) calloc( (field.number)*(field.number), sizeof(double) );
    if(fr1 == NULL){fprintf(stderr,"Allocation error, fr1 in lens_forvard\n");
			exit(1);}

      fi1 = (double *) calloc( (field.number)*(field.number), sizeof(double) );
    if (field.imaginary == NULL){fprintf(stderr,"Allocation error, fi1 in lens_forvard\n");
			exit(1);}
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
	  int i_i, j_i;
	  long ik1;
	  i_i=field.number-i+1;
	  j_i=field.number-j+1;
	  ik1 = (i_i-1)*field.number + j_i -1;
	  ik=(i-1)*field.number +j-1;
            fr1[ik] = field.real[ik1]/ampl_scale;
            fi1[ik] =field.imaginary[ik1]/ampl_scale; 

        }

    }
for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
	 
	  ik=(i-1)*field.number +j-1;
            field.real[ik]= fr1[ik];
            field.imaginary[ik] = fi1[ik];

        }
}
free(fi1);
free(fr1);
    }
    /*
            fprintf(stderr,"%e %e %e %e %e\n",f1,f,  z1, field.size,1./field.double1);
            */
    write_field();


}

#include "fft_prop.c"




void error_print(char *arr)
{
    fprintf(stderr,"\n%s propagates the field to \
distance Z [units you use]\nin VARIABLE COORDINATE  system",arr);
    fprintf(stderr,"\n\nUSAGE:  ");
    fprintf(stderr,"%s F Z, where  F is \
the focal lenght of the input lens,\n\
Z is the distance to propagate\n\n",arr);

}







