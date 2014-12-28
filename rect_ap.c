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

    int i,j,i2;
    double sx,sy,dx,x,y, x0,y0,x_shift, y_shift, angle,cc,ss;
    long ik1;



    /* Processing the command line argument  */

    if (argc<2 || argc >6){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&sy);
    sx=sy;
    x_shift=y_shift=angle=0.;
    if(argc == 3) sscanf(argv[2],"%le",&sx);
    if(argc == 4) {
        sscanf(argv[2],"%le",&sx);
        sscanf(argv[3],"%le",&y_shift);
    }
    if(argc == 5){
        sscanf(argv[2],"%le",&sx);
        sscanf(argv[3],"%le",&y_shift); 
        sscanf(argv[4],"%le",&x_shift);
    }
     if(argc == 6){
        sscanf(argv[2],"%le",&sx);
        sscanf(argv[3],"%le",&y_shift); 
        sscanf(argv[4],"%le",&x_shift);
	sscanf(argv[5],"%le",&angle);
	angle *= -Pi2/360.;
    }
    read_field();

    dx =field.size/(field.number);
    i2=field.number/2+1;

    /* Cuttitng the aperture      */

    for (i=1;i<=field.number ;i++)
        for (j=1;j<=field.number ;j++){

            ik1=(i-1)*field.number+j-1;
            x0=(i-i2)*dx-x_shift;
            y0=(j-i2)*dx-y_shift;

	    cc=cos(angle);
	    ss=sin(angle);
	    x=x0*cc+y0*ss;
	    y=-x0*ss+y0*cc;  
            if(fabs(x) > sx/2. || fabs(y) > sy/2. ) {
                field.real[ik1]=0.;
                field.imaginary[ik1]=0.;
            }

        }
    write_field();


}


void error_print(char *arr) { 

    fprintf(stderr,"\n%s filters the field \
through a square aperture \n",arr);


    fprintf(stderr,"\nUSAGE: %s X [Y DX DY Angle], where X and Y are \n\
the sides of the  aperture in [units you use],\n\
DX and DY are the shifts of the aperure center\n\
Angle is the rotation of the aperture in degrees of arc\n\
The aperture is first shifted, then roatated\n\n",arr);

}







