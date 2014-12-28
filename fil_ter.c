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
FILE *fr;
long ik;
double *sum, cc, smax, ssum, fi, cab,sab,phi,ccc;
int i,j, bin_ind, ind;
char * first_string, c_int;
unsigned char b_in;
ind=0;

 if (argc< 4 || argc >5 ){
        error_print(argv[0]);
        exit(1);
    }

read_field();


first_string = (char *) calloc (110,1);

if((fr=fopen(argv[3],"r"))==NULL){
fprintf(stderr,"error opening file %s, exiting \n", argv[3]);
exit(1);} 
bin_ind =0;

i=0;
  while((c_int=getc(fr)) != '\n') {
    first_string[i]=c_int;
    i++;
    if(i>100) i=100;
}
if((strstr(first_string, "P2")) != NULL ) {
/*  fprintf(stderr,"Portable ASCII graymap detected \n"); */
}
else if((strstr(first_string, "P5")) != NULL ){
/*  fprintf(stderr,"Portable binary graymap detected \n"); */
 bin_ind=1;
}    
else goto l1;
  
do{i=0;
  while((c_int=getc(fr)) != '\n') {
    first_string[i]=c_int;
    i++;
    if(i>100) i=100;}
/*  fprintf(stderr, "%s\n", first_string); */
  } while (first_string[0] == '#');
  
do{ i=0;
  while((c_int=getc(fr)) != '\n') {
    first_string[i]=c_int;
    i++;
    if(i>100) i=100;}
/*  fprintf(stderr, "%s\n", first_string); */
  } while (first_string[0] == '#');
goto l2;
  l1:
rewind(fr);
  l2:
  if ((strstr(argv[1], "in"))!= NULL) {
  if ((strstr(argv[2], "su"))!=NULL){
/*  int subst here */
sum  = (double *) calloc( (field.number)*(field.number), sizeof(double) );
if(sum == NULL){fprintf(stderr,"fil_ter int subst: Allocation error, exiting\n");
exit(1);}

   smax=0;
   ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
	  if(bin_ind) {
	   if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
	  fprintf(stderr,"Error reading portable bitmap\n");
	  exit(1);}
	   sum[ik]= (double) b_in;
	  }
	  else{
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter int subst: end of input file reached, exiting\n");
	      exit(1);}
	      sum[ik]=fi;
	    	  }
    
	  
if(sum[ik] < 0 && ind == 0){ 
	      fprintf(stderr,"fil_ter int subs  warning: the\
 intensity is negative\n"); ind=1;
	      }
	    if(smax < sum[ik]) smax=sum[ik];
	    ik++;
	}
    }
    
    if (argc != 5) smax=1.;

    ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
	    ssum=sqrt(sum[ik]/smax);
	    phi=phase(field.imaginary[ik], field.real[ik]); 
	    field.real[ik]=ssum*cos(phi);
	    field.imaginary[ik] = ssum*sin(phi);
	    ik++;
	  }
      }
  free(sum);
  
}
  else { 
/* int mult here */
sum  = (double *) calloc( (field.number)*(field.number), sizeof(double) );
if(sum == NULL){fprintf(stderr,"fil_ter int subst: Allocation error, exiting\n");
exit(1);}

   smax=0;
   ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){

 if(bin_ind) {
	   if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
	  fprintf(stderr,"Error reading portable bitmap\n");
	  exit(1);}
	   sum[ik]= (double) b_in;
	  }
	  else{
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter int subst: end of input file reached, exiting\n");
	      exit(1);}
	      sum[ik]=fi;
	  } 
             

    
	  
if(sum[ik] < 0 && ind == 0){
	      fprintf(stderr,"fil_ter int mult  warning: the\
 intensity is negative\n"); ind = 1;
	      }
	    if(smax < sum[ik]) smax=sum[ik];
	    ik++;
	}
    }

if (argc != 5) smax=1.;

    ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
	    ssum=sqrt(sum[ik]/smax);
	     
	    field.real[ik] *= ssum;
	    field.imaginary[ik] *= ssum;
	    ik++;
	  }
      }
  free(sum);
  }
}
else { 
if (strstr(argv[2], "su")!=NULL){
/*  pha  subst here */
ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){
 
    if(bin_ind) {
	   if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
	  fprintf(stderr,"Error reading portable bitmap\n");
	  exit(1);}
	   fi = (double) b_in;
	  }
	  else{
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter pha subst: end of input file reached, exiting\n");
	      exit(1);}
	      
	  } 
/*          if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter pha subst: end of input file reached, exiting\n");
	      exit(1);}*/

	    cab=cos(fi);
            sab=sin(fi);
	    ccc=sqrt(field.real[ik]*field.real[ik]\
+field.imaginary[ik]*field.imaginary[ik]);
field.real[ik]= ccc*cab;
field.imaginary[ik] = ccc*sab;
ik ++;
	}
}
	       
}
  else { /* pha  mult here */
  
    ik=0;
    for (i=1;i<= field.number;i += 1){
        for (j=1;j <= field.number;j += 1){ 
              if(bin_ind) {
	   if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
	  fprintf(stderr,"Error reading portable bitmap\n");
	  exit(1);}
	   fi = (double) b_in;
	  }
	  else{
              if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter int subst: end of input file reached, exiting\n");
	      exit(1);}
	      
	  } 
/*  if ((fscanf(fr,"%le",&fi))==EOF){
	      fprintf(stderr,"fil_ter pha mult: end of input file reached, exiting\n");
	      exit(1);}*/

	    cab=cos(fi);
            sab=sin(fi);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
	    ik++;

	  }
      }
  }

}


fclose(fr);
write_field();
free(first_string);
}

void error_print(char *arr){

fprintf(stderr,"\n%s allows for operations on the phase/intensity\n",arr);

fprintf(stderr,"\nUSAGE: %s C1 C2 F [N], where\n\
C1 is character constant with valid values int and pha\n\
C2 is character constant with valid values mult and subst\n\
F is the name of a file containing intensity/phase mask\n\
It may be portable graymap or anymap (*.pgm, *.pnm) with grayscale data.\n",arr);
fprintf(stderr,"\
The number of values in F must be the same as the grid dimension\n\
N is (any) optional argument indicating that the intensity mask should\n\
be normalized before applying\n\n\
Examples:\n\
%s int mult aa: filter THROUGH  intensity mask from file aa\n\
%s pha subst aa: substitute the phase with phase taken from file aa\n\
%s pha mult aa: filter the field through phase filter aa\n\
%s int subst aa: substitute intensity with one taken from file aa\n\
%s int subst aa haha: substitute intensity with a normalized\n\
one taken from file aa\n\n",arr,arr,arr,arr,arr);
}

