/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    


#include "pipes.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h> 
void main(int argc, char *argv[]){
    void error_print();
    int i,j;
 unsigned int seed;
    double ampl;

    
    long ik;
    /* Processing the command line argument  */
    if (argc< 2 || argc >4 ){
        error_print(argv[0]);
        exit(1);
    }
    sscanf(argv[2],"%le",&ampl);
   seed=(unsigned int) time((unsigned int) 0);
/*    fprintf(stderr,"random seed = %d \n", seed); */
    if(argc == 4) sscanf(argv[3],"%d",&seed);
    srand(seed);
    
       read_field();
    if ((strstr(argv[1], "ph"))!= NULL){

    ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){ 
            double fi, cab,sab, my_rnd, cc;
	    my_rnd= ((double) rand()) / ((double) RAND_MAX)-0.5;
	    fi=my_rnd*ampl;
	    cab=cos(fi);
            sab=sin(fi);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
	    ik++;

	  }
      }
    }
    else if ((strstr(argv[1], "in"))!= NULL){
      double maxint =0.;
       ik=0;
       for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){ 
            double   cc;
	    cc=field.imaginary[ik]*field.imaginary[ik] + \
field.real[ik]*field.real[ik];
	    if (cc > maxint) maxint=cc; 
	    ik++;
	  }
      }
      ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){ 
            double  cab, sab, my_rnd, cc, phi;
	    my_rnd= ((double) rand()) / ((double) RAND_MAX);
	    phi=phase(field.imaginary[ik],field.real[ik]);
	     cc=(field.imaginary[ik]*field.imaginary[ik] + \
field.real[ik]*field.real[ik])+ maxint*ampl*my_rnd;
	    cc=sqrt(cc);
	    cab=cos(phi);
	    sab=sin(phi);
            field.imaginary[ik]=cc*cab;
            field.real[ik]=cc*sab;
	    ik++;
	  }
      }
    }
    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s  filters the field through \
a random  mask\n",arr);

    fprintf(stderr,"\nUSAGE: %s I A [S]\n\
where I may be <int> for the intensity or <pha> for the phase\n\
A is the amplitude of phase noise in radians or;\n\
A is the amplitude of the intensity noise,\n\
normalized to the maximum intensity\n\
S (optional) is an unsigned int seed, if omitted S\n\
is taken from the computer clock. \n\
Warning: the filter changes the integral intensity \n\n",arr);


}

