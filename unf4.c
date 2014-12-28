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
double  *x , dum,   koeff;
int  i ,  nn,j,  ii, ii1, ii2, ii3, ubbf3();
long ik;

 if (argc<2 || argc >2){
        error_print(argv[0]);
        exit(1);
    }
 sscanf(argv[1],"%le",&koeff);
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
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1+1)*nn +j-1;
    ii3=(i-1)*nn +j-1+1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
} 
for(i=nn; i>=2; i--){
  for (j=1; j<=nn-1; j++){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1-1)*nn +j-1;
    ii3=(i-1)*nn +j-1+1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
} 
for(i=1; i<=nn-1; i++){
  for (j=nn; j>=2; j--){
/*    fprintf(stderr, "%d %d \n",i,j);*/
    ii1=(i-1)*nn+j-1;
    ii2=(i-1+1)*nn +j-1;
    ii3=(i-1)*nn +j-1-1;
    ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
      }
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





int ubbf3( a1, a2, a3)
double *a1, *a2, *a3; { 
  double s[10], aa2, aa3;
  int ic,  ik;
 
double d5(), Pi2;
int i_min();
Pi2=Pi*2.;

  l1:
  
ik=1;

  for(aa2= *a2 - Pi2; aa2<= *a2 + Pi2+0.1; aa2 += Pi2){
    for(aa3= *a3 - Pi2; aa3<= *a3 + Pi2+0.1; aa3 += Pi2){
      s[ik]=d5( *a1, aa2, aa3 );
/*      fprintf(stderr,"%ld %e\n", ik, s[ik]);*/
      ik++;
    }
}
 
ic=i_min(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9]);
/*fprintf(stderr,"%d %e %e %e %e %e\n", ic, *a5, *a1, *a2, *a3, *a4); */

if (ic==1){
  *a2 -= Pi2;
  *a3 -= Pi2;
  goto l1;
}
if (ic==2){
   *a2 -= Pi2;
  goto l1;
}
if (ic==3){
   *a2 -= Pi2;
   *a3 += Pi2;
  goto l1;
}
if (ic==4){
   *a3 -= Pi2;
  goto l1;
}
if (ic==5){
   goto l2;
}
if (ic==6){
   *a3 += Pi2;
  goto l1;
}
if (ic==7){
   *a2 += Pi2;
   *a3 -= Pi2;
  goto l1;
}
if (ic==8){
   *a2 += Pi2;
    goto l1;
}
if (ic==9){
   *a2 += Pi2;
   *a3 += Pi2;
    goto l1;
}


  l2: return ic;
}


double d5(a1,a2,a3)
double a1,a2,a3;
  {
return fabs(a1-a2)+fabs(a1-a3);
}



int i_min(a1,a2,a3, a4, a5, a6, a7, a8, a9)
double a1,a2,a3,a4,a5,a6,a7,a8,a9;
  { int ii;
     double dmin;
if (a2 < a1) {dmin=a2; ii=2;}
else {dmin=a1;         ii=1;}
if (a3<dmin) {dmin=a3; ii=3;}
if (a4<dmin) {dmin=a4; ii=4;}
if (a5<dmin) {dmin=a5; ii=5;}
if (a6<dmin) {dmin=a6; ii=6;}
if (a7<dmin) {dmin=a7; ii=7;}
if (a8<dmin) {dmin=a8; ii=8;}
if (a9<dmin) {dmin=a9; ii=9;}  
  return ii;
}




void error_print(char *arr) { 

    fprintf(stderr,"\n%s unfolds the phase saved in a Gnuplot format\n",arr);


    fprintf(stderr,"\nUSAGE: %s  K < phase_in > phase_out \n\
where K is coefficient to multiply the phase, may be\n\
necessary to translate the phase into dimensional units.\n\n",arr); 

}









