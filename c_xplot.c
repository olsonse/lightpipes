#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define Pi 3.141592654


double aa[4096];
void error_print();

void main(int argc, char *argv[]){
    double *x, dum, xmax= -1e23, xmin=1e23;
    int i ,  nn,  j;
    long ik, n_vert, n_pol, ikk;
    FILE *fr, *f_out ;

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
            fprintf(stderr,"f_xplot: Allocation error while reading, exiting\n");
            exit(1);
        }
        x[ik]=dum;	 
	/* fprintf(stderr," %e \n", x[ik]); */
        ik++;
    }
    nn=(int) sqrt( (double) ik);
    n_vert=nn*nn;
    n_pol=0;
   fclose(fr); 

   ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	  if(xmin > x[ik]) xmin=x[ik];
	  if(xmax < x[ik]) xmax=x[ik];
	 ik++;
	   }
    } 

xmax -= xmin;

f_out=fopen(argv[2],"w");
   ik=0;
    ikk=0;
    for(i=1; i<=nn-1; i++){
        for (j=1; j<=nn-1; j++){
	  double v1, v2, v3,v4;
	  v1=x[ik];
	  v2=x[ik+1];  
	  v3=x[ik+nn];
	  v4=x[ik+nn+1];
	    if(v1 != 0 && v2 != 0 && v3 != 0 && v4 != 0  ){
	    	     ikk++;
	    }
 
	 ik++;
	   }
    } 

fprintf(f_out,"%ld %ld ",ik,ikk);
    ik=0;
    for(i=1; i<=nn; i++){
        for (j=1; j<=nn; j++){
	  fprintf(f_out,"%g %g %g \n",(double) i, (double) j, x[ik]*nn/xmax);
	 ik++;
	   }
    } 
    n_vert=ik;
 
    ik=0;
    ikk=0;
    for(i=1; i<=nn-1; i++){
        for (j=1; j<=nn-1; j++){
	  double v1, v2, v3,v4;
	  v1=x[ik];
	  v2=x[ik+1];  
	  v3=x[ik+nn];
	  v4=x[ik+nn+1];
	    if(v1 != 0 && v2 != 0 && v3 != 0 && v4 != 0  ){
	     fprintf(f_out,"4 \t %ld %ld %ld %ld \n",ik, ik+1, ik+1+nn, ik+nn);
	     ikk++;
	    }
 
	 ik++;
	   }
    } 

    fclose(f_out);

    free(x);

    fclose(fr);
}

void error_print(char *arr) { 

    fprintf(stderr,"\n%s converts the intensity/phase file into xplot file \n", arr);

    fprintf(stderr,"\nUSAGE: %s F  M_F  where F is the file in gnuplot format,\n\
M_F is the xplot readable .off file\n\n", arr);

}







