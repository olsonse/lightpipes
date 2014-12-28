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
double  *x , dum,   koeff, rad;
int  i ,  nn,j,  ii, ii1, ii2,   ubbf2(), ilong=2, icount;
long ik ;

 if (argc<2 || argc >3){
        error_print(argv[0]);
        exit(1);
    }
 sscanf(argv[1],"%le",&koeff);
sscanf(argv[2],"%le",&rad);
x  = (double *) calloc( 1, sizeof(double) );
ik=0;

while ((scanf("%le", &dum))!= EOF){

x = (double *) realloc ( x, sizeof(double)*(ik+1));
if(x == NULL){fprintf(stderr,"unfold: Allocation error while reading, exiting\n");
exit(1);}
x[ik]=dum;
ik++;
}

nn=(int) sqrt( (double) ik);
for(i=1; i<=nn-1; i++){
  for (j=1; j<=nn-1; j++){
    int ii, ij;
    ii = i-nn/2+1;
    ij=j-nn/2+1; 
    ii1=(i-1)*nn+j-1;
    if (ii*ii+ij*ij > rad*rad*nn*nn/4) x[ii1]=0.;
}
}

for (icount=1; icount <= ilong; icount ++){

for(i=nn/2; i<=nn-1; i++){
  for (j=nn/2; j<=nn-1; j++){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1)*nn +j-1+1;
    ii=ubbf2(&x[ii1],&x[ii2]);
      }
} 

for(i=nn/2; i>=2; i--){
  for (j=nn/2; j<=nn-1; j++){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1-1)*nn +j-1;
    
    ii=ubbf2(&x[ii1],&x[ii2]);
      }
} 
for(i=nn/2; i<=nn-1; i++){
  for (j=nn/2; j>=2; j--){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1+1)*nn +j-1;
    
    ii=ubbf2(&x[ii1],&x[ii2]);
      }
} 
for(i=nn/2; i>=2; i--){
  for (j=nn/2; j>=2; j--){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1)*nn +j-1-1;
    
    ii=ubbf2(&x[ii1],&x[ii2]);
      }
} 
/*
for(i=nn; i>=2; i--){
  for (j=1; j<=nn-1; j++){

    ii1=(i-1)*nn+j-1;
    ii2=(i-1-1)*nn +j-1;
    ii3=(i-1)*nn +j-1+1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
} 

for(i=nn; i>=2; i--){
  for (j=nn; j>=2; j--){

    ii1=(i-1)*nn+j-1;
    ii2=(i-1-1)*nn +j-1;
    ii3=(i-1)*nn +j-1-1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
} 

for(i=1; i<=nn-1; i++){
  for (j=nn; j>=2; j--){
    ii1=(i-1)*nn+j-1;
    ii2=(i-1+1)*nn +j-1;
    ii3=(i-1)*nn +j-1-1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
} 
*/
}


for(i=1; i<=nn; i++){
  for (j=1; j<=nn; j++){
    ik=(i-1)*nn+j-1;
    printf("%e\n", x[ik]*koeff);
	   }
printf("\n");
}
free(x);

}







int ubbf2( a1, a2)
double *a1, *a2; { 
double s[5], aa2;  
int ic,  ik;
 
double d3(), Pi2;
int i_min();
Pi2=Pi*2.;

l1:
  
ik=1;

  for(aa2= *a2 - Pi2; aa2<= *a2 + Pi2+0.1; aa2 += Pi2){
 
      s[ik]=d3( *a1, aa2 );
/*      fprintf(stderr,"%ld %e\n", ik, s[ik]);*/
      ik++;
    }

 
ic=i_min(s[1],s[2],s[3]);

if (ic==1){
  *a2 -= Pi2;
  goto l1;
}

if (ic==2){
  goto l2;
}
if (ic==3){
   *a2 += Pi2;
  goto l1;
}



  l2: return ic;
}

double d3(a1,a2)
double a1,a2;
  {
return fabs(a1-a2);
}



int i_min(a1,a2,a3)
double a1,a2,a3;
  { int ii;
     double dmin;
if (a2 < a1) {dmin=a2; ii=2;}
else {dmin=a1;         ii=1;}
if (a3<dmin) {dmin=a3; ii=3;}
  return ii;
}


void error_print(char *arr) { 

    fprintf(stderr,"\n%s unfolds the phase saved in a Gnuplot format\n",arr);


    fprintf(stderr,"\nUSAGE: %s  K rad < phase_in > phase_out \n\
where K is coefficient to multiply the phase, may be\n\
necessary to translate the phase into dimensional units.\n\
rad is the radius of blackening aperture. relative to grid size\n",arr); 

}










