/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    


#include "pipes.h"



void main(int argc, char *argv[]){

    void error_print();
    int i,j,imax, istep, ii, jj;
    
    FILE *fr;
    double /* im,im0,*/ al, im, re;
    long ik1;
    /* Processing the command line argument  */
    if (argc<2 || argc >4){
	error_print(argv[0]);
	exit(1);
    }

    read_field();


    imax=64;
    

    if (argc >= 3){

    if((strstr(argv[2], "sam"))!= NULL ) {imax=field.number;}
    else{
    if((sscanf(argv[2],"%d",&imax))!=0){}
    else imax=field.number;
    }

    }
    al=0.01;
    if(argc ==4) {
      sscanf(argv[3],"%le",&al);}

    istep=1;
    if (imax>field.number) imax=field.number;
    if(field.number/imax >1) {istep= field.number/imax;
    imax=field.number/istep;}



/*    im=0;
    for (i=1;i<=field.number;i+=istep){
	for (j=1;j<=field.number;j+=istep){
	    ik1=(i-1)*field.number+j-1;
	    im0=field.real[ik1]*field.real[ik1]+\
field.imaginary[ik1]*field.imaginary[ik1];
	    if (im < im0) im=im0;
	}
    }
*/

#ifdef _DJGPP_
 setmode(fileno(stdout), O_BINARY);

  fr=fopen(argv[1],"wb");
#else
	fr=fopen(argv[1],"w");
#endif  


    /* writing the phase    */
    if(istep != 1){
    for (i=1;i < field.number;i+=istep){
	for (j=1;j< field.number;j+=istep){
	  re=im=0.;
	  for(ii=i; ii<i+istep; ii++){
	    for(jj=j; jj<j+istep; jj++){
	      
	    ik1=(ii-1)*field.number+jj-1;
	    re += field.real[ik1];
	    im += field.imaginary[ik1];
	    }
	  }
	
	    
	  fprintf(fr,"%e\n",phase(im,re));
	  

	}
	fprintf(fr,"\n");

    }
    }
    else { 
      for (i=1;i<= field.number;i+=istep){
	for (j=1;j<= field.number;j+=istep){
	  re=im=0.;
	  ik1=(i-1)*field.number+j-1;
	  re = field.real[ik1];
	  im = field.imaginary[ik1];
	  fprintf(fr,"%e\n",phase(im,re));
	}
	fprintf(fr,"\n");
    }
    }
    fclose(fr);

    write_field();

}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s writes phase of the field  into file F\n\n",arr);

    fprintf(stderr,"USAGE: %s F [N L], where F is the filename,\n\n\
optional parameter N is the number of points in the grid, \n\
N=64 if omitted in command line,\n\
N equals to grid sampling if you pass the word \"same\"\n\
L is the level of intensity at which the phase is put out as zero\n\
default L=0.01\n\n",arr);
fprintf(stderr,"The phase can be plotted with gnuplot splot command\n\n");



}










