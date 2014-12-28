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
    double split, asplit, bsplit, csplit;
    int i,j;

    FILE *fr;
    long ik1;
    /* Processing the command line argument  */
#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);
#endif


    if (argc< 2 || argc >3 ){
	error_print(argv[0]);
	exit(1);
    }
    if (argc == 3) sscanf(argv[2],"%le",&split);
    else split=0.5;
    if(split == 0.) split=0.000001;



    asplit=sqrt(split);
    bsplit=sqrt(1-split);
    csplit=bsplit/asplit;


    read_field();

    for (i=1;i<=field.number ;i++)
	for (j=1;j<=field.number ;j++){
	    ik1=(i-1)*field.number+j-1;

	    field.real[ik1] *= asplit;
	    field.imaginary[ik1] *=asplit;
	}

    { 
	long size;

#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);
  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif

	if(fwrite ((char *)&(field.number), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"b_split Error while writing FIELD.number\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.size), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"b_split Error while writing FIELD.size\n");
	    exit(1);
	}
	if(fwrite ((char *)&(field.lambda), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"b_split Error while writing FIELD.lambda\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.int1), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.int1\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.int2), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.int2\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.int3), sizeof(int), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.int3\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.double1), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.double1\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.double2), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.double2\n");
	    exit(1);
	}

	if(fwrite ((char *)&(field.double3), sizeof(double), 1, fr) != 1){
	    fprintf(stderr,"Error while writing FIELD.double3\n");
	    exit(1);
	}


	size=(field.number);
	size *=size;
	size=size*sizeof(double);
	if(fwrite (field.real, (long unsigned) size, 1, fr)!=1){
	    fprintf(stderr,"b_split Error while writing FIELD.real\n");
	    exit(1);
	}
	if(fwrite (field.imaginary, (long unsigned) size, 1, fr)!=1){
	    fprintf(stderr,"b_split Error while writing FIELD.imaginary\n");
	    exit(1);
	}
	fclose(fr);
    }


    for (i=1;i<=field.number ;i++)
	for (j=1;j<=field.number ;j++){
	    ik1=(i-1)*field.number+j-1;

	    field.real[ik1] *= csplit;
	    field.imaginary[ik1] *= csplit;
	}


    write_field();
}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s  splits the beam, producing two: \
\n", arr);

    fprintf(stderr,"\nUSAGE: %s F [A], where F is the output \
filename,\
\noptional parameter A is the splitting coefficient, A=0.5 \
if omitted\n\n", arr);

    /*fprintf(stderr,"\n\n");*/


}







