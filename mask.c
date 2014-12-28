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
    double *x, *x1, dum;
    double  lev;
    int i ,  nn, nn1, j;
    long ik;
    FILE *fr, *fr1;

    if (argc< 2 || argc >4 ){
        error_print(argv[0]);
        exit(1);
    }

sscanf(argv[3],"%le",&lev);

    fr=fopen(argv[1],"r");
    x = (double *) calloc( 1, sizeof(double) );
    ik=0;
    while ((fscanf(fr,"%le", &dum))!= EOF){
        x = (double *) realloc ( x, sizeof(double)*(ik+1));
        if(x == NULL){
            fprintf(stderr,"mask: Allocation error while reading, exiting\n");
            exit(1);
        }
        x[ik]=dum;	 
	/* fprintf(stderr," %e \n", x[ik]); */
        ik++;
    }
    nn=(int) sqrt( (double) ik);
 
   x1 = (double *) calloc( 1, sizeof(double) );
 
        fr1=fopen(argv[2],"r");

        ik=0;

        while ((fscanf(fr1,"%le", &dum))!= EOF){

            x1 = (double *) realloc ( x1, sizeof(double)*(ik+1));
            if(x1 == NULL){
                fprintf(stderr,"mask: Allocation error while reading, exiting\n");
                exit(1);
            }
            x1[ik]=dum;	 
	    /* fprintf(stderr," %e \n", x1[ik]); */
            ik++;
        }
        nn1=(int) sqrt( (double) ik);
	if (nn1 != nn)  fprintf(stderr,"mask: dimensions of files are different, exiting\n");
	fclose(fr1);
    

/* extracting one function from the another */
    ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	 if(x1[ik] > lev) printf("%e\n",x[ik]);
	 else printf("%e\n",0.);
	 ik++;
	   }
	printf("\n");
    } 

    free(x);
    free(x1);
    fclose(fr);
}

void error_print(char *arr) { 

    fprintf(stderr,"\n%s prints to the output the first file masked by the second file\n",arr);

    fprintf(stderr,"\nUSAGE: %s  P1 P2 level, where P1 and P1 are \n\
names of datafiles containing the phase distributions\n\
level is the masking level, the value from the first file will be put out\n\
if the correspondent value of the second file is greater than level\n\n",arr);

}







