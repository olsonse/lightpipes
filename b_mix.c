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
    double f_real, f_imag, fieldsize, fieldlambda, fielddouble1;
    int i,j, fint2;

    FILE *fr;
    long ik1;
    int fieldnumber;

    /* Processing the command line argument  */
    if (argc< 2 || argc >2 ){
	error_print(argv[0]);
	exit(1);
    }



    read_field();
    fint2=field.int2;


    { 
#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"rb");
#else
	fr=fopen(argv[1],"r");
#endif
       
	if(fread (&fieldnumber, sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.number\n");
	    exit(1);
	}

	if(fieldnumber != field.number){
	    fprintf(stderr,"b_mix: You can not mix grids with \
different dimensions!\n");
	    exit(1);
	}

	if(fread (&(fieldsize), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.size\n");
	    exit(1);
	}

	if(fieldsize != field.size){
	    fprintf(stderr,"b_mix: You can not mix grids with \
different sizes!\n");
	    exit(1);
	}

	if(fread (&(fieldlambda), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.lambda\n");
	    exit(1);
	}

	if(fieldlambda != field.lambda){
	    fprintf(stderr,"b_mix: You can not mix grids with \
different wavelengths!\n");
	    exit(1);
	}


	if(fread (&(field.int1), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.int1\n");
	    exit(1);
	}

  
	if(fread (&(field.int2), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.int2\n");
	    exit(1);
	}
        if(fint2 >field.int2) field.int2=fint2;

	if(fread (&(field.int3), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.int3\n");
	    exit(1);
	}

	if(fread (&(fielddouble1), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.double1\n");
	    exit(1);
	}

	if(fielddouble1 != field.double1){
	    fprintf(stderr,"b_mix: You can not mix grids originated \
from different coordinate systems!\n");
	    exit(1);
	}

	if(fread (&(field.double2), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.double2\n");
	    exit(1);
	}

	if(fread (&(field.double3), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while reading FIELD.double3\n");
	    exit(1);
	}



	for (i=1;i<=field.number ;i++)
	    for (j=1;j<=field.number ;j++){
		ik1=(i-1)*field.number+j-1;
		if(fread (&f_real, sizeof(double), 1, fr)!=1){
		    fprintf(stderr,"Error while reading FIELD.real in b_mix \n");
		    exit(1);
		}



		field.real[ik1] += f_real;

	    }
	for (i=1;i<=field.number ;i++)
	    for (j=1;j<=field.number ;j++){ 
		ik1=(i-1)*field.number+j-1;
		if(fread (&f_imag, sizeof(double), 1, fr)!=1){
		    fprintf(stderr,"Error while reading FIELD.imaginary in b_mix \n");
		    exit(1);
		}



		field.imaginary[ik1] += f_imag;
	    }

	fclose(fr);
    }





    write_field();
}




void error_print( char *arr)
{
    fprintf(stderr,"\n%s mixes two coherent fields:  \
\n",arr);

    fprintf(stderr,"\nUSAGE: %s  F < F1, where F and F1 are the input \
filenames\
\n\n", arr);

    /*fprintf(stderr,"\n\n");*/

}

