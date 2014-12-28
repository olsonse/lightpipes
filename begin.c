/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

/* This filter forms the initial grid filled with
   unity field                               
*/

#include "pipes.h"


void main( int argc, char *argv[]){

void errorprint();
    int n_grid,i,j;
    long ik;
    double size_grid, lambda;





    /* Processing the command line argument  */

    if (argc < 3 || argc >4) {
	errorprint(argv[0]);
    }


    /* reading the data from a command line */


    sscanf(argv[1],"%le",&size_grid);
    sscanf(argv[2],"%le",&lambda);
    n_grid=256;
    if(argc>3) sscanf(argv[3],"%d",&n_grid);

   

    newfield(n_grid, size_grid, lambda,0,0,0,0.,0.,0.);


    /*  Here the initial field is formed   */
    ik=0;
    for (i=1;i<=field.number ;i++)
	for (j=1;j<=field.number ;j++)

	{
	    field.real[ik]=1.;
	    field.imaginary[ik]=0.;
	    ik++;
	}


    write_field(); 
    /*    return (0);
			*/
}

void errorprint(char *arr){


    fprintf(stderr,"\n%s forms the initial data structure for use \
by all following programs\n", arr);
    fprintf(stderr,"\nUSAGE: begin B C [A],  where:                       \n\n");
    fprintf(stderr,"A is the grid dimension, grid of AxA \
points will be formed (if omitted, A=256) \n");
    fprintf(stderr,"B is the grid size, the real size of the grid will\
 be put as BxB [units you use]\n");
    fprintf(stderr,"C is the wavelength in [units you use]\n\n");
    exit(1);

}

























