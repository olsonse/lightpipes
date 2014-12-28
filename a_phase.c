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
    double *x, *x1, dum, sig, sig1,di,dj;
    double f_f, f_i, f_j, c1, ci,cj,cii,cij,cjj, tilti, tiltj, shift;
    int i ,  nn, nn1, j;
    long ik;
    FILE *fr, *fr1;

    if (argc< 2 || argc >3 ){
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
    x1 = (double *) calloc( 1, sizeof(double) );
    if(argc == 2){
      x1 = (double *) realloc( x1, nn*nn* sizeof(double) );
      ik=0;
      for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	  x1[ik]=0;
	  ik++;
	}
      }
    }	  

/* the second file */
     else if (argc == 3) {
        fr1=fopen(argv[2],"r");

        ik=0;

        while ((fscanf(fr1,"%le", &dum))!= EOF){

            x1 = (double *) realloc ( x1, sizeof(double)*(ik+1));
            if(x1 == NULL){
                fprintf(stderr,"a_phase: Allocation error while reading, exiting\n");
                exit(1);
            }
            x1[ik]=dum;	 
	    /* fprintf(stderr," %e \n", x1[ik]); */
            ik++;
        }
        nn1=(int) sqrt( (double) ik);
	if (nn1 != nn)  fprintf(stderr,"a_phase: dimensions of files are different, exiting\n");
	fclose(fr1);
    }

/* extracting one function from the another */
    ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	 dum = x[ik]-x1[ik];
	  x1[ik] =dum;
	 ik++;
	   }
    } 
/* calculation of  the shift and tilts x */
    f_f=f_i=f_j=c1=ci=cj=cii=cij=cjj=0.;
    shift=tilti=tiltj=0.;

    ik=0;
    for(i=1; i<=nn; i++){
      di=(double) i;
        for (j=1; j<=nn; j++){
	  dj=(double) j;

	  f_f += x[ik];
	  f_i += di*x[ik];
	  f_j += dj*x[ik];
	  c1  += 1.;
	  ci  += di;
	  cj  += dj;
	  cii += di*di;
	  cjj += dj*dj;
	  cij += di*dj;
	  ik++;
	   }
    } 

shift = -(cjj*ci*f_i-cjj*f_f*cii-ci*cij*f_j-cj*cij*f_i+cj*f_j*cii+f_f*cij*cij)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

tilti = (-ci*cjj*f_f-cj*cj*f_i+cj*cij*f_f+cj*ci*f_j+cjj*c1*f_i-cij*c1*f_j)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

tiltj = -(c1*cij*f_i-c1*f_j*cii-ci*cj*f_i-ci*cij*f_f+ci*ci*f_j+f_f
*cj*cii)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

fprintf(stderr," %e %e %e\n", shift, tilti, tiltj);

/* extraction of tilts x */
 ik=0;
    for(i=1; i<=nn; i++){ 
      di=(double) i;
      for (j=1; j<=nn; j++){
	 dj=(double) j;
	x[ik] = x[ik] - di*tilti - dj*tiltj - shift;
	 ik++;
	   }
    } 

/* calculation of  the shift and tilts x1 */
    f_f=f_i=f_j=c1=ci=cj=cii=cij=cjj=0.;
    shift=tilti=tiltj=0.;

    ik=0;
    for(i=1; i<=nn; i++){
      di=(double) i;
        for (j=1; j<=nn; j++){
	  dj=(double) j;

	  f_f += x1[ik];
	  f_i += di*x1[ik];
	  f_j += dj*x1[ik];
	  c1  += 1.;
	  ci  += di;
	  cj  += dj;
	  cii += di*di;
	  cjj += dj*dj;
	  cij += di*dj;
	  ik++;
	   }
    } 

shift = -(cjj*ci*f_i-cjj*f_f*cii-ci*cij*f_j-cj*cij*f_i+cj*f_j*cii+f_f*cij*cij)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

tilti = (-ci*cjj*f_f-cj*cj*f_i+cj*cij*f_f+cj*ci*f_j+cjj*c1*f_i-cij*c1*f_j)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

tiltj = -(c1*cij*f_i-c1*f_j*cii-ci*cj*f_i-ci*cij*f_f+ci*ci*f_j+f_f
*cj*cii)/(-cj*cj*cii+cjj*c1*cii-ci*ci*cjj+2*cj*cij*ci-cij*cij*c1);

fprintf(stderr," %e %e %e\n", shift, tilti, tiltj);

/* extraction of tilts x1 */
 ik=0;
    for(i=1; i<=nn; i++){ 
      di=(double) i;
      for (j=1; j<=nn; j++){
	 dj=(double) j;
	x1[ik] = x1[ik] - di*tilti - dj*tiltj - shift;
	 ik++;
	   }
    } 
/* dispersions */


    sig=sig1=0;
    ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	  sig += x[ik]*x[ik];
	  sig1 += x1[ik]*x1[ik];
	  ik ++;
	}
    }
    sig /= (double) ik;
    sig1 /= (double) ik;



    fprintf(stderr," Sigma square before %e  after %e ratio %e\n", sig, sig1, sig1/sig);
   fprintf(stderr," Sigma  before %e  after %e ratio %e\n", sqrt(sig), sqrt(sig1), sqrt(sig1/sig));

 ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	  printf ("%e\n",x1[ik]);
	  ik ++;
	}
	printf("\n");
    }

    free(x);
    free(x1);
    fclose(fr);
}

void error_print(char *arr) { 

    fprintf(stderr,"\n%s prints to the output the info about the phase distribution\n",arr);

    fprintf(stderr,"\nUSAGE: %s  P1 [P2], where P1 and P1 are \n\
names of datafiles containing the phase distributions\n\
the program will output to the stdout the phase map with shift extracted\n\
if two files are given the program will extract the second phase\n\
from the first one and then output the difference to the stdout\n\n", arr);

}







