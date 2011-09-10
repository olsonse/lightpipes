/*--------------------------------------------------------------*/
/*      (C) modifications Gleb Vdovin 1993-1996                 */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
/*--------------------------------*-C-*---------------------------------*
 * \file fftn.h.
 * ---------------------------------------------------------------------*
 * Re[]:	real value array
 * Im[]:	imaginary value array
 * nTotal:	total number of complex values
 * nPass:	number of elements involved in this pass of transform
 * nSpan:	nspan/nPass = number of bytes to increment pointer
 *		in Re[] and Im[]
 * isign:	exponent: +1 = forward  -1 = reverse
 * scaling:	normalizing constant by which the final result is *divided*
 *	scaling == -1, normalize by total dimension of the transform
 *	scaling <  -1, normalize by the square-root of the total dimension
 *
 * ----------------------------------------------------------------------
 * See the comments in the code for correct usage!
 */

#ifndef _FFTN_H
#define _FFTN_H

/** Free-up allocated temporary storage after finished all the Fourier
 * transforms.
 */
extern void fft_free (void);

/** Compute fft with double precision.
 *
 * @param ndim
 *      Total number dimensions.
 * @param dims
 *      Vector of array sizes.
 *      If NDIM is zero then DIMS must be zero-terminated
 *
 * @param Re
 *      Holds the real component of the data, and returns the resulting real
 *      Fourier coefficient.  Multidimensional data *must* be allocated
 *      contiguously.  There is no limit on the number of dimensions. 
 * @param Im
 *      Holds the imaginary component of the data, and returns the resulting
 *      imaginary Fourier coefficient.  Multidimensional data *must* be
 *      allocated contiguously.  There is no limit on the number of
 *      dimensions.  
 *
 * @param isign
 *      Sign of the complex exponential (ie, forward or inverse FFT)
 *	the magnitude of ISIGN (normally 1) is used to determine the
 *	correct indexing increment (see below).
 *
 * @param scaling
 *      Normalizing constant by which the final result is *divided*
 *	if scaling == -1, normalize by total dimension of the transform
 *	if scaling <  -1, normalize by the square-root of the total dimension
 *
 * example:<br>
 * tri-variate transform with Re[n1][n2][n3], Im[n1][n2][n3]
 *
 *	int dims[3] = {n1,n2,n3}
 *	fftn (3, dims, Re, Im, 1, scaling);
 *
*/
extern int fftn (int ndim, const int dims[], double Re[], double Im[],
		 int isign, double scaling);

/* * Compute fft with float precision. */
extern int fftnf (int ndim, const int dims[], float Re[], float Im[],
		  int isign, double scaling);
#endif	/* _FFTN_H */
