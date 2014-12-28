/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define Pi 3.141592654


double aa[4096];
void error_print();

void main(int argc, char *argv[]){
double  *x , dum,   koeff, amax;
int  i ,  nn,j;
long ik ;

 if (argc<2 || argc >2){
        error_print(argv[0]);
        exit(1);
    }
 sscanf(argv[1],"%le",&koeff);

x  = (double *) calloc( 1, sizeof(double) );
ik=0;

while ((scanf("%le", &dum))!= EOF){

x = (double *) realloc ( x, sizeof(double)*(ik+1));
if(x == NULL){fprintf(stderr,"negate: Allocation error while reading, exiting\n");
exit(1);}
x[ik]=dum;
ik++;
}

nn=(int) sqrt( (double) ik);

amax=0;
for(i=1; i<=nn; i++){
  for (j=1; j<=nn; j++){
    ik=(i-1)*nn+j-1;
    if ( x[ik] > amax) amax=x[ik];

      }
} 

for(i=1; i<=nn; i++){
  for (j=1; j<=nn; j++){
    ik=(i-1)*nn+j-1;
    printf("%e\n", koeff*(amax-x[ik]));
	   }
printf("\n");
}
free(x);

}




void error_print(char *arr) { 

    fprintf(stderr,"\n%s reads intensity distribution and outputs the same negative\n",arr);


    fprintf(stderr,"\nUSAGE: %s  K < int_in > int_out \n\
where K is coefficient to multiply the negated intensity, may be\n\
necessary to scale the distribution.\n\n",arr); 

}










