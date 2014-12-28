/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

void forvard(zz)
double zz;
{ 
    int i,j, ii, ij, n12;
    long ik, ir;
    double z,z1,cc;
    double sw, sw1, bus, abus, pi2, cab, sab;
    pi2=2.*3.141592654;
    z=fabs(zz);
 
    ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]*ii*ij;
            field.imaginary[ik] =field.imaginary[ik]*ii*ij; 
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }

    if(zz>=0.) fft3(1);
    else fft3(-1);

    /*
        Spatial filter, (c)  Gleb Vdovin  1986:  
                                    */ 
    if(zz >= 0.){
    z1=z*field.lambda/2.;
    n12=field.number/2;
    ik=0;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){ 
            sw=((i-n12-1)/field.size);
            sw *= sw;
            sw1=((j-n12-1)/field.size);
            sw1 *= sw1;
            sw += sw1; 
            bus=z1*sw;
            ir = (long) bus;
            /*	    if (ir > 1e8) fprintf(stderr,"%ld\n",ir);  */
            abus=pi2*(ir- bus);

            cab=cos(abus);
            sab=sin(abus);
            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
            field.real[ik]=cc;
            ik++;
        }
    }
}
    else { 
      z1=z*field.lambda/2.;
      n12=field.number/2;
      ik=0;
      for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){ 
            sw=((i-n12-1)/field.size);
            sw *= sw;
            sw1=((j-n12-1)/field.size);
            sw1 *= sw1;
            sw += sw1; 
            bus=z1*sw;
            ir = (long) bus;
            /*	    if (ir > 1e8) fprintf(stderr,"%ld\n",ir);  */
            abus=pi2*(ir- bus);

            cab=cos(abus);
            sab=sin(abus);
            cc=field.real[ik]*cab + field.imaginary[ik]*sab;
            field.imaginary[ik]= field.imaginary[ik]*cab-field.real[ik]*sab;
            field.real[ik]=cc;
            ik++;
        }
    }
    
}

 
  
    if(zz >= 0.) fft3(-1);
    else fft3(1);   
    
	ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]*ii*ij;
            field.imaginary[ik] =field.imaginary[ik]*ii*ij;
     
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }
    
}

 
void fft3(ind)
int ind;
{int dims[2];
dims[0]=dims[1]=field.number;
/* fprintf(stderr,"%d %d \n", dims[0], dims[1]); */
fftn(2, dims, field.real, field.imaginary, ind, (double) field.number);
}



















