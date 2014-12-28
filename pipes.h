/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 

#ifdef _DJGPP_
#include <fcntl.h>
#endif



#define Pi 3.14159265358979323844
#define Pi2 6.28318530717958647688



#ifdef _DJGPP_
extern void setmode(int,int);
#endif


void newfield();
void write_field();
void read_field();



/* The structure FIELD contains the characteristics
of the light beam: number of points along the side 
of a  square grid, wavelength and side length of the
square grid, then two huge arrays of Re and Im data 
*/

typedef struct{ 
    int number;
    double size, lambda;
    int int1,int2,int3;
    double double1,double2,double3;
    double *real, *imaginary; 
}
FIELD;
FIELD field;
double phase();
char  * pass_string, p_pass[10]; 

void read_field()
/* program read_field reads the FIELD from the standard input */
{

   
    unsigned long  size;

 



    
#ifdef _DJGPP_
  setmode(fileno(stdin), O_BINARY);
#endif

 if(fread (&(field.number), sizeof(int), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.number\n" );
	  exit(1);}

    if( field.number/2 != (float) (field.number)/2.){
	  fprintf(stderr,"Sorry, number of points must be even, stopping\n");
	  exit(1);}



 if(fread (&(field.size), sizeof(double), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.size\n");
	  exit(1);}
 if(fread (&(field.lambda), sizeof(double), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.lambda\n");
	  exit(1);}

 if(fread (&(field.int1), sizeof(int), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.int1\n");
	  exit(1);}
 if(fread (&(field.int2), sizeof(int), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.int2\n");
	  exit(1);}





 if(fread (&(field.int3), sizeof(int), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.int3\n");
	  exit(1);}



 if(fread (&(field.double1), sizeof(double), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.double1\n");
	  exit(1);}
 if(fread (&(field.double2), sizeof(double), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.double2\n");
	  exit(1);}
 if(fread (&(field.double3), sizeof(double), 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD.double3\n");
	  exit(1);}



    newfield(field.number, field.size, field.lambda, field.int1, \
field.int2, field.int3, field.double1, field.double2, field.double3);
    size=(field.number);
    size *=size;
    size=size*sizeof(double);
	if(fread (field.real, size, 1, stdin) != 1){
	  fprintf(stderr,"Error while reading FIELD real\n");
	  exit(1);}
 
	if(fread (field.imaginary, size, 1, stdin)!=1){
     fprintf(stderr,"Error while reading FIELD imaginary\n");
	 exit(1);}

}


void write_field()
/* program write_field writes field to the standard output */
{

    unsigned long  size;

#ifdef _DJGPP_
  setmode(fileno(stdout), O_BINARY);
#endif    

 if(fwrite ((char *)&(field.number), sizeof(int), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.number\n");
	  exit(1);}

 if(fwrite ((char *)&(field.size), sizeof(double), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.size\n");
	  exit(1);}
if(fwrite ((char *)&(field.lambda), sizeof(double), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.lambda\n");
	  exit(1);}

 if(fwrite ((char *)&(field.int1), sizeof(int), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.int1\n");
	  exit(1);}

 if(fwrite ((char *)&(field.int2), sizeof(int), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.int2\n");
	  exit(1);}

 if(fwrite ((char *)&(field.int3), sizeof(int), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.int3\n");
	  exit(1);}

 if(fwrite ((char *)&(field.double1), sizeof(double), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.double1\n");
	  exit(1);}

 if(fwrite ((char *)&(field.double2), sizeof(double), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.double2\n");
	  exit(1);}

 if(fwrite ((char *)&(field.double3), sizeof(double), 1, stdout) != 1){
	  fprintf(stderr,"Error while writing FIELD.double3\n");
	  exit(1);}

	size=field.number;
	size *= size;
	size *= sizeof(double);
    
	if(fwrite ((char *) field.real, size, 1, stdout)!=1){
	  fprintf(stderr,"error writing field.real\n");
	    exit(1);}



       if( fwrite ((char *) field.imaginary, size, 1, stdout)!=1){
	  fprintf(stderr,"error writing field.imaginary\n");
	    exit(1);}
       

   free (field.real);
   free (field.imaginary);
    
}

/* newfield allocates memory */

void newfield(number,size,lambda, int1,int2,int3, doub1, doub2, doub3)
unsigned int number;
double size, lambda, doub1, doub2, doub3;
int int1, int2, int3;
 

{

    field.number=number;
    field.size=size;
    field.lambda=lambda;
    field.int1=int1;
    field.int2=int2;
    field.int3=int3;
    field.double1=doub1;
    field.double2=doub2;
    field.double3=doub3;

    field.real = (double *) calloc( (number)*(number), sizeof(double) );
    if(field.real == NULL){fprintf(stderr,"Allocation error, Real\n");
			exit(1);}

    field.imaginary = (double *) calloc( (number)*(number), sizeof(double) );
    if (field.imaginary == NULL){fprintf(stderr,"Allocation error, Imaginary\n");
			exit(1);}
}


double phase(y,x)
double y,x;
{    double pp=0.; 
     if(x==0.){
       if (y>0) pp= 0.5*Pi;
       if (y==0)pp=0.;
       if (y<0) pp= -0.5*Pi;
     }
     else {
if(y !=0) pp= atan2 (y,x);
else pp=0.;
     }
return pp;

 
}














