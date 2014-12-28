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
    int ind, ii,ij, iiij, i, j;
    long ik;
    

    /* Processing the command line argument  */
    if (argc!=2){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%d",&ind);
    
    
 

    read_field();
     field.int1 += ind;
    if( field.int1 != 0 ){
      
     ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
	  iiij=ii*ij;
            field.real[ik] *= iiij;
            field.imaginary[ik] *= iiij;
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }
   }

    
    fft3(ind);

 if(field.int1 == 0){  
 
     ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
	  iiij=ii*ij;
            field.real[ik] *= iiij;
            field.imaginary[ik] *= iiij;
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }
   }    
    write_field();


}






void error_print(char *arr)
{
    fprintf(stderr,"\n%s  performs FFT of the field distribution\n",arr);
    fprintf(stderr,"\n\nUSAGE:  ");
    fprintf(stderr,"%s IND , where IND is direction,\n\
IND=1 for forward, -1 for inverse \n\n",arr);

}

void fft3(ind)
int ind;
{int dims[2];
dims[0]=dims[1]=field.number;
/* fprintf(stderr,"%d %d \n", dims[0], dims[1]); */
fftn(2, dims, field.real, field.imaginary, ind, (double) field.number);
}
