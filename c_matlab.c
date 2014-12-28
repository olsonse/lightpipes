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
    double *x, dum;
    int i ,  nn,  j;
    long ik;
    FILE *fr ;

    if (argc< 2 || argc >2 ){
        error_print(argv[0]);
        exit(1);
    }



    fr=fopen(argv[1],"r");
    x = (double *) calloc( 1, sizeof(double) );
    ik=0;
    while ((fscanf(fr,"%le", &dum))!= EOF){
        x = (double *) realloc ( x, sizeof(double)*(ik+1));
        if(x == NULL){
            fprintf(stderr,"unfold: Allocation error while reading, exiting\n");
            exit(1);
        }
        x[ik]=dum;	 
	/* fprintf(stderr," %e \n", x[ik]); */
        ik++;
    }
    nn=(int) sqrt( (double) ik);
   fclose(fr); 

/* otput matlab */
printf ("\nA=[");

    ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	 if (j != nn ) printf("%g ",x[ik]);
	 else{
	   if ( i != nn ) printf("%g\n",x[ik]);
	   else printf("%g",x[ik]);
	 }
	 ik++;
	   }
    } 
printf("];\n");

    free(x);

    fclose(fr);
}

void error_print(char *arr) { 

    fprintf(stderr,"\n%s converts the intensity/phase file into matlab m-file \n", arr);

    fprintf(stderr,"\nUSAGE: %s F > M_F  where F is the file in gnuplot format, M_F is the matlab m-file\n\n", arr);

}







