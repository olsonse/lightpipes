/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    


#include "pipes.h"




void main(int argc, char *argv[]){
    double * int1;
    double int2;
    double * phase1;
    double phase2;
    
    void error_print();
    int i,j, jj;
    double dx;
    FILE *fr;
    long ik1;

    /* Processing the command line argument  */

    if (argc< 2 || argc >2 ){
	error_print(argv[0]);
	exit(1);
    }

    read_field();


    dx=field.size/(field.number-1.);

    /* allocating memory for output arrays  */

    int1=(double *) calloc(field.number, sizeof(double));
    if(int1 == NULL){
	fprintf(stderr,"Allocation error in cros_out , int1\n");
	exit(1);
    }

    phase1=(double *) calloc(field.number, sizeof(double));
    if(int1 == NULL){
	fprintf(stderr,"Allocation error in cros_out , int1\n");
	exit(1);
    }



    /* writing the intensity into arrays  */

    i=field.number/2+1;
    for (j=1;j<=field.number;j += 1){ 
	ik1=(i-1)*field.number+j-1;
	jj=j-1; 
	int1[jj]=field.real[ik1]*field.real[ik1]+field.imaginary[ik1] *field.imaginary[ik1]; 
	phase1[jj]=phase(field.imaginary[ik1],field.real[ik1]);

    }
#ifdef _DJGPP_
 setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif


    j=field.number/2+1;
    for (i=1;i<=field.number;i += 1){
	double cc;
	ik1=(i-1)*field.number+j-1;
	jj=i-1;
	cc=dx*(i-field.number/2-1);

	int2=field.real[ik1]*field.real[ik1]+field.imaginary[ik1] *field.imaginary[ik1]; 
	phase2=phase(field.imaginary[ik1], field.real[ik1]);
	fprintf(fr," %e %e %e %e %e\n", cc, int1[jj], int2, phase1[jj], phase2);

    }


    fclose(fr);


    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s  writes X and Y cross sections of intensity \
and phase distributions into file F\n",arr);

    fprintf(stderr,"USAGE: %s  F, where F is the output filename,\n\n\
first column   :    coordinate\n\n\
second column  :    intensity distribution along X coordinate,\n\
section through the grid center \n\
\nthird column   :    intensity distribution along Y coordinate,\n\
section through the grid center\n\
\nfourth column  :    phase distribution along X coordinate, \n\
section through the grid centern\n\
\nfifth  column  :    phase distribution along X coordinate, \n\
section through the grid center\n\
    \n",arr);
    fprintf(stderr,"The distributions  can be plotted with gnuplot plot command\n\n");


}


