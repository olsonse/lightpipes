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
    int i,j, imax, istep;
	double dx,dx2;
    FILE *fr;
    long ik1;
    /* Processing the command line argument  */
    if (argc< 2 || argc >3 ){
	error_print(argv[0]);
	exit(1);
    }
   
    read_field();
    imax=64;

    if (argc == 3){

    if((strstr(argv[2], "sam"))!= NULL ) {imax=field.number;}
    else{
    if((sscanf(argv[2],"%d",&imax))!=0){}
    else imax=field.number;
    }

    }
	dx=field.size/(field.number-1.);
	dx2=dx*dx;

    istep=1;
    if(imax>field.number)imax=field.number;
    if(field.number/imax >1) {istep= field.number/imax;
    imax=field.number/istep;}


#ifdef _DJGPP_
 setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif
   


    /* writing the intensity     */

    for (i=1;i<= field.number;i += istep){
	for (j=1;j <= field.number;j += istep){ 
	    double sum;



	    sum=0;
/*	    for (ii=i; ii<i+istep; ii++)
		for(jj=j; jj<j+istep; jj++){ */
		    ik1=(i-1)*field.number +j- 1; 
		    sum += field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
/*		}
	    sum=sum/(istep*istep);
	    */
	    fprintf(fr,"%e\n", sum );

	}
	fprintf(fr,"\n");

    }
    fclose(fr);


    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s writes intensity \
distribution into file F\n",arr);

    fprintf(stderr,"\nUSAGE: %s F [N], where F is the output filename,\n\
optional parameter N is the number of points in the plotting grid, \n\
N=64 if omitted in command line\n\
N equals to grid sampling if you pass the word \"same\" \n",arr);
    fprintf(stderr,"The intensity can be plotted with gnuplot splot command\n\n");


}


