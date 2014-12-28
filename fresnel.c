/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <math.h> 
#include "pipes.h" 

void main(int argc, char *argv[]){
double z;
void fresnel();
void error_print();
    
    
    /* Processing the command line argument  */

    if (argc != 2 ){
        error_print(argv[0]);
        exit(1);
    }

    /* reading the data from a command line */

sscanf(argv[1],"%le",&z);
read_field();    
    
    fresnel(z);
 /*   fprintf(stderr,"%d %e %e \n",field.number, field.lambda, field.size); */
    write_field();
 
}

void error_print(char *arr)
{
    fprintf(stderr,"\n%s propagates the field to \
distance Z [units you use]\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s  Z , where  Z is the distance to propagate\n\
A convolution method of calculation of the \n\
diffraction integral is implemented here.\n\
There is no reflection an the grid borders.\n\
On the same grid the function uses 8 times more RAM then forvard.\n\n",arr);

}

#include "frsn_prop.c"




















