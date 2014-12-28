/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/

/* changes the contrast constant gamma for gnuplot formatted file */

    
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void main(int argc, char *argv[]){
void error_print();
int n, i, j;
long ik;
double gamma, *fi, amax,cc;
   if (argc< 2 || argc >3 ){
        error_print(argv[0]);
        exit(1);
    }
 if (argc >= 2) sscanf(argv[1],"%d",&n);
    else n=256; 
 if (argc >= 3) sscanf(argv[2],"%le",&gamma);
    else gamma=1;

fi  = (double *) calloc( n*n , sizeof(double) );
if(fi == NULL){fprintf(stderr,"Allocation error, exiting\n");
exit(1);}

ik=0;
amax=0;
for (i=1; i<=n; i += 1){ 
for (j=1; j <= n; j += 1){ 
              if ((scanf("%le",&cc))==EOF) {
	      fprintf(stderr,"end of input file reached, exiting\n");
	      exit(1);}
	      fi[ik]=cc;
	     if(amax < fi[ik]) amax=fi[ik];
	      ik++;
	}
    }

 ik=0;
for (i=1; i<=n; i++){ 
for (j=1;j <= n;j += 1){ 
  printf("%e \n", amax * pow((fi[ik]/amax),1./(gamma+0.00001)) );
  ik++;
}
printf("\n");
}

}

void error_print (char *arr) { 

    fprintf(stderr,"\n%s changes the contrast constant gamma for gnuplot formatted file \n", arr);

    fprintf(stderr,"\nUSAGE: %s n G < file1 > file2\n\
where n is grid sampling, G is the gamma value\n\n", arr);

}






