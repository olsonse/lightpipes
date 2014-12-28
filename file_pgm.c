/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include "pipes.h"
#include <math.h>
#define GAMMA 2.0

void main(int argc, char *argv[]){

    void error_print();
    int i,j,ii,jj,imax, istep,i0, i_i, max_val;

    FILE *fr;
    long ik1;
    float max_int=0, gamma;
    /* Processing the command line argument  */
    if (argc< 2 || argc >5 ){
	error_print(argv[0]);
	exit(1);
    }
    read_field();
  imax=field.number;
  if (argc >= 3){

    if((strstr(argv[2], "sam"))!= NULL ) {imax=field.number;}
    else{
    if((sscanf(argv[2],"%d",&imax))!=0){}
    else imax=field.number;
    }
  }



    if (argc >= 4) sscanf(argv[3],"%e",&gamma);
    else gamma=GAMMA;
    if (argc >= 5) sscanf(argv[4],"%d",&max_val);
    else max_val=255;





  /*    if (((double) field.number)/((double) imax) != field.number/imax) imax=field.number;
   */
   

    istep=1;
    if(imax>field.number)imax=field.number;
   if(field.number/imax >1) {
istep= field.number/imax;
imax=(int)ceil ((float)field.number/(float) istep);
    }

#ifdef _DJGPP_
 setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif  
  
/*    fprintf(stderr,"%d\n", istep); */
    if(istep != 1){
    /* header of the PNM  file */
    fprintf(fr,"P2\n");
    fprintf(fr,"#Creator: LightPipes (C) 1993-1996, Gleb Vdovin\n");
    fprintf(fr,"%d %d\n", imax-1, imax-1);
    fprintf(fr,"%d\n", max_val);
    }
    else{ /* header of the PNM  file */
    fprintf(fr,"P2\n");
    fprintf(fr,"#Creator: LightPipes (C) 1993-1996, Gleb Vdovin\n");
    fprintf(fr,"%d %d\n", imax, imax);
    fprintf(fr,"%d\n", max_val);
    }

    if( istep != 1){

    for (i=1 ; i<= field.number-istep; i +=  istep){
	for (j=1;j <= field.number-istep;j += istep){
	    double sum;
	    sum=0;
	    for (ii=i; ii<=i+istep; ii++)
		for(jj=j; jj<=j+istep; jj++){
		    ik1=(ii-1)*field.number +jj- 1;
		    sum += field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
		}
	    sum=sum/(istep*istep);
	    if(sum>max_int) max_int=sum;
	}
    }


    i_i=1;
    for (i=1; i<= field.number-istep; i += istep){
       for (j=1;j <= field.number-istep;j += istep){ 
	    double sum;
	    sum=0;
	    for (ii=i; ii<=i+istep; ii++)
		for(jj=j; jj<=j+istep; jj++){
		    ik1=(ii-1)*field.number +jj- 1;
		    sum +=field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
		}
	    sum=sum/(istep*istep);
	    i0=(int) floor(pow((sum/max_int),1./(gamma+0.0001))*max_val);
	   /* i0=(int)  (sum/max_int)*255;*/
	    fprintf(fr,"%d ", i0);
	    i_i++;
	    if (i_i == 40){
	      fprintf(fr,"\n");
		i_i=1;
	    }

	}


    }
    }
    else{
      
      for (i=1; i<= field.number; i++){
       for (j=1 ;j <= field.number; j++ ){
	double sum;
	ik1=(i-1)*field.number +j- 1;
	 sum =field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
	if(sum>max_int) max_int=sum;
       }
      }
      i_i=1;
      for (i=1; i<= field.number; i++ ){
       for (j=1; j <= field.number; j++ ){
	double sum;
	ik1=(i-1)*field.number +j- 1;
	 sum =field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
	i0=(int) floor(pow((sum/max_int),1./(gamma+0.00001))*max_val);
	   
	    fprintf(fr,"%d ", i0);
	    i_i++;
	    if (i_i == 40){
	      fprintf(fr,"\n");
		i_i=1;
	    }
       }
      }
    }


    fclose(fr);


    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s writes intensity distribution\
into *.pgm file F\n",arr);

    fprintf(stderr,"\nUSAGE: %s F [N, G, MAX], where F is the output filename,\n\
optional parameter N is the grid size, N=128 if omitted in command line\n\
N equals to grid sampling if you pass the word \"same\"\n\
G is the Gamma parameter, [0.1...10], higher G gives better\n\
contrast in low intensities, default G=2.0\n\
MAX is the number of gray levels, default MAX=255\n\n",arr);
    fprintf(stderr,"Output file F can be processed with netpbm package \n\n");


}

#undef GAMMA






















