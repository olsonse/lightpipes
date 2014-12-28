/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

 
#include "pipes.h"
#include <math.h>
#include <string.h>

void main(int argc, char *argv[]){
    void error_print();

    int i,j, n2;
    double sum, sum1r, sum1i, sum1, sum2, dx,dx2, s_p, x,y, x_c, y_c;
    double sum1x, sum1y;
    long ik1;



    /* Processing the command line argument  */

    if (argc!= 2){
        error_print(argv[0]);
        exit(1);
    }


    read_field();

    dx =field.size/(field.number);
    dx2 = dx*dx;
    n2=field.number/2+1;



    /* Calculating the power */
    sum=sum1r=sum1i=sum2=0.;
    ik1=0;
    for (i=1;i<=field.number ;i++){
        for (j=1;j<=field.number ;j++){



            
	      s_p=(field.real[ik1]*field.real[ik1]+ \
field.imaginary[ik1]*field.imaginary[ik1]);
	    sum2 += s_p;
	    sum += sqrt(s_p);
            sum1r += field.real[ik1];
            sum1i += field.imaginary[ik1];
            ik1++;
        }
    }
    sum1=(sum1r*sum1r+sum1i*sum1i);


    if (sum == 0) {

fprintf(stderr,"Strehl: Zero beam power, program terminated\n");
        exit(1);
			      }

if(strstr(argv[1], "y")!= NULL)fprintf(stderr,"Strehl: ratio= %e energy= %e\n",sum1/sum/sum, sum2*dx2);

/* Calculating the center of gravity: */
 sum=sum1r=sum1i=sum2=0.;
 ik1=0;
    for (i=1;i<=field.number ;i++){
      y=(i-n2)*dx;
        for (j=1;j<=field.number ;j++){
	  x=(j-n2)*dx;
	  sum2=(field.real[ik1]*field.real[ik1]\
+field.imaginary[ik1]*field.imaginary[ik1]);
	  sum1r += sum2*x;
	  sum1i += sum2*y;
	  sum += sum2;

            ik1++;
        }
    }

    x_c=sum1r/sum;
    y_c=sum1i/sum;

   fprintf(stderr,"Center_of_gravity: x= %e y= %e\n", x_c, y_c);
  

/* Calculating moments of the distribution */
 sum1r=sum1x=sum1y=0.;
 ik1=0;
    for (i=1;i<=field.number ;i++){
      double y_y_c;
      y=(i-n2)*dx;
      y_y_c=y-y_c;
      
        for (j=1;j<=field.number ;j++){
	  double temp_int, x_x_c;
	  x=(j-n2)*dx;
	  x_x_c=x-x_c;
	  temp_int = (field.real[ik1]*field.real[ik1]\
+field.imaginary[ik1]*field.imaginary[ik1]);
	  sum1r += temp_int*(x_x_c*x_x_c+y_y_c*y_y_c);
	  sum1x += temp_int*(x_x_c*x_x_c);
	  sum1y += temp_int*(y_y_c*y_y_c);

            ik1++;
        }
    }


   fprintf(stderr,"Standard deviation:  S_r=%e S_x= %e S_y= %e\n", sqrt(sum1r/sum), sqrt(sum1x/sum), sqrt(sum1y/sum));
 fprintf(stderr,"Grid size: %e, Grid sampling: %d\n", field.size, field.number); 


    write_field();


}


void error_print(char *arr) { 

    fprintf(stderr,"\n%s: prints the general info to the stderr\n",arr);


    fprintf(stderr,"\n%s y, y prevents arrival of this message\n\n",arr);




}







