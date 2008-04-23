/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*      (C) C++ port and modifications Spencer Olson 2006       */
/*      For bugs in C++ port, send bug reports to               */
/*      olsonse at umich.edu                                    */
/*                                                              */
/*--------------------------------------------------------------*/
    

#include <lightpipes/Field.h>
#include <math.h>

void error_print(const char *arr = NULL);

int main(int argc, const char *argv[]) {

    double K;

    /* Processing the command line argument  */

    if (argc<2 || argc >2){
        error_print(argv[0]);
        exit(1);
    }
    /* reading the data from a command line */
    sscanf(argv[1],"%le",&K);
    double kk=sqrt(K);

    Field * field = Field::read (  );
    (*field) *= kk;
    field->write (  );
    delete field;

    return 0;
}


void error_print(const char *arr) {
    fprintf(stderr,"\n%s filters the field through \
intensity attenuator K\n", arr);

    fprintf(stderr,"\nUSAGE: %s K, where K \
is the coefficient \nof attenuation or amplification (must be positive)\n\n", arr);


}


