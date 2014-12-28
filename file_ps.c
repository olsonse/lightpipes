/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    


#include "pipes.h"

#define GAMMA 2.0

void main(int argc, char *argv[]){

    void error_print();
    int i,j,ii,jj,imax, istep,i0;

    FILE *fr;
    long ik1;
    float max_int=0, gamma;
    /* Processing the command line argument  */
    if (argc< 2 || argc >4 ){
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





  /*    if (((double) field.number)/((double) imax) != field.number/imax) imax=field.number;
   */
    

    istep=1;
    if(imax>field.number)imax=field.number;
    if(field.number/imax >1) 
{istep= field.number/imax;
imax=(int)ceil ((float)field.number/(float) istep);}

    /* fprintf(stderr,"%s: istep= %d\n", argv[0], istep); */
#ifdef _DJGPP_
 setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif  
  


    /* header of the postscript file */
    fprintf(fr,"%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(fr,"%%%%BoundingBox: 0 0 256 256  \n");
    fprintf(fr,"%%%%Creator: LightPipes (C) 1993 Gleb Vdovin, beta version 1995 \n");
 /*   fprintf(fr,"%%%%Pages: 1\n");*/
    fprintf(fr,"%%%%EndComments\n");
    fprintf(fr,"%%%%EndProlog\n\n");
/*    fprintf(fr,"%%%%Page: 1 1\n\n");*/

    fprintf(fr,"/origstate save def\n20 dict begin\n");
    fprintf(fr,"/picstr %d string def\n", (imax-1));
    fprintf(fr,"256  256   scale\n");
    fprintf(fr,"%d %d  8\n", (imax-1), (imax-1));
    fprintf(fr,"[ %d 0 0 %d 0 %d]\n",(imax-1), -(imax-1), (imax-1));
    fprintf(fr,"{currentfile\npicstr readhexstring pop}\nimage\n ");
    

    for (i=1; i<= field.number-istep; i +=  istep){
	for (j=1; j <= field.number-istep; j += istep){
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



    for (i = 1; i<= field.number-istep; i +=  istep){
       for (j=1;j <= field.number-istep;j += istep){ 
	    double sum;
	    sum=0;
	    for (ii=i; ii<=i+istep; ii++)
		for(jj=j; jj<=j+istep; jj++){
		    ik1=(ii-1)*field.number +jj- 1;
		    sum +=field.real[ik1] *field.real[ik1]+ field.imaginary[ik1] *field.imaginary[ik1];
		}
	    sum=sum/(istep*istep);
	    i0=(int) floor(pow((sum/max_int),1./(gamma+0.00001))*255);
	   /* i0=(int)  (sum/max_int)*255;*/
	    fprintf(fr,"%02x", i0);

	}


    }

   fprintf(fr,"\n\nshowpage\n%%%%Trailer\n");

    fclose(fr);


    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s writes intensity \
into PostScript file F\n",arr);

    fprintf(stderr,"\nUSAGE: %s F [N, G], where F is the output filename,\n\
optional parameter N is the grid size, N=128 if omitted in command line\n\
N equals to grid sampling if you pass the word \"same\"\n\
G is the Gamma parameter, [0.1...10], higher G gives better\n\
contrast in low intensities, default G=2.0\n\n",arr);
    fprintf(stderr,"Output file F can be printed with any postscript device\n\n");


}

#undef GAMMA






















