/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
#include "pipes.h"
double Zernike();
void Zer();


void main(int argc, char *argv[]){
    void error_print();
    double R,A;
    int n,m;

    /* Processing the command line argument  */
    if (argc!=5){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%d",&n);
    sscanf(argv[2],"%d",&m);
    sscanf(argv[3],"%le",&R);
    sscanf(argv[4],"%le",&A);
    /*fprintf(stderr," %d %d %e %e\n", n, m, R, A); */
    read_field();

    Zer(n,m,R,A); 

    write_field();

}

void Zer(n,m,R,A)
int n,m;
double R,A;
{ 
    int i,j,n2;
    long ik;
    double x, y, dx, fi, cab, sab, cc, cc1,rho,phi;

    n2=field.number/2;

    dx=field.size/field.number;

    ik=0;
    /*     fprintf(stderr,"%le %le\n", tx,ty);*/ 
    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx;
        for (j=1;j<=field.number; j++){
            y=(j-n2-1)*dx;
            rho=sqrt((x*x+y*y)/(R*R));

            phi=phase(y,x);


            fi= A*Zernike(n,m,rho,phi);
            cab=cos(fi);
            sab=sin(fi);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            cc1=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
	    field.imaginary[ik]=cc1;
            ik++;
        }

    }
}




void error_print(char *arr)
{
    fprintf(stderr,"\n%s  introduces arbitrary Zernike aberration into the field distribution\n",arr);
    fprintf(stderr,"\nUSAGE:  ");
    fprintf(stderr,"%s  n m R A, where  n and m  are the integer orders - \n\
see Born and Volf p. 465, sixth (corrected) edition, Pergamon, 1993,\n\
R is the radius at which the phase amplitude reaches A (in radians)\n\n",arr);

}


/***************************************************************/
/* Zernike polynomial 

   +-m
  R    as in Born and Volf, p. 465, sixth edition, Pergamon
     n

The implementation have not been optimized for speed.
 
*/


double Zernike(n,m,rho,phi)
int n,m;
double rho, phi;
{
    int s, int_sign, mm, ncheck, ind;
    double sum, product, factorial();
    if(n<0){
        fprintf(stderr,"Zernike: n must be >0; |m| must be less or equal than n\n\
if n is odd then m must be odd,\n\
if n is even then m must be even\n");
        exit(1);
    }
    ind=0;
    for(ncheck=n; ncheck>=-n; ncheck -= 2){
        if (ncheck == m) ind=1;
    }
    if(ind == 0){
        fprintf(stderr,"Zernike: n must be >0; |m| must be less or equal than n\n\
if n is odd then m must be odd,\n\
if n is even then m must be even\n");
        exit(1);
    }

    mm=(int) fabs(m);
    sum=0;
    int_sign=1;
    for (s=0; s<= (int)((n-mm)/2); s++){
        if(n-2*s != 0) product=pow( rho, (double)(n-2*s) );
        else product =1;
        product *= factorial(n-s)*int_sign;
        product /= factorial(s)*factorial(((n+mm)/2)-s)*factorial(((n-mm)/2)-s);
        sum += product;
        int_sign = -int_sign;
    }

    if(m<=0) return sum*cos(m*phi);
    else return sum*sin(m*phi); 
}



/* Factorial function */
double factorial(n)
int n;
{
    double product;

    if (n<0) {
        fprintf(stderr,"factorial: argument is negative, exiting \n");
        exit(1);
    }
    if (n==0) return 1.;
    else{ 
        product =1;
        while(n>=1){
            product *= n;
            --n;
        }
        return product;
    }
}




/*****************************************************************/
/*****************        END of Zernike     *********************/ 





