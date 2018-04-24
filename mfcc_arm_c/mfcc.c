/*
 * Based on https://github.com/jsawruk/libmfcc and https://www.nayuki.io/res/how-to-implement-the-discrete-fourier-transform/dft.c and https://batchloaf.wordpress.com/2013/12/07/simple-dft-in-c/
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libmfcc.h"

#define N 16
#define PI2 6.2832

void compute_dft_from_signal(const double signal[], double power_spectrum[]) {
    // time and frequency domain data arrays
    int n, k;             // indices for time and frequency domains
    float Xre[N], Xim[N]; // DFT of x (real and imaginary parts)
     
    // Calculate DFT of x using brute force
    for (k=0 ; k<N ; ++k)
    {
        // Real part of X[k]
        Xre[k] = 0;
        for (n=0 ; n<N ; ++n) Xre[k] += signal[n] * cos(n * k * PI2 / N);
         
        // Imaginary part of X[k]
        Xim[k] = 0;
        for (n=0 ; n<N ; ++n) Xim[k] -= signal[n] * sin(n * k * PI2 / N);
         
        // Power at kth frequency bin
        power_spectrum[k] = Xre[k]*Xre[k] + Xim[k]*Xim[k];
    }
}

int main(void) {
	double spectrum[8192];

	// Pointer to the sample data file
	FILE *sampleFile;

	// Index counter - used to keep track of which data point is being read in
	int i = 0;

	// Determine which MFCC coefficient to compute
	unsigned int coeff;

	// Holds the value of the computed coefficient
	double mfcc_result;

	// Initialize the spectrum
	memset(&spectrum, 0, sizeof(spectrum)); 
	
	// Open the sample spectrum data	
	sampleFile = fopen("sample.dat","rb");
	
	// Read in the contents of the sample file
	while(fscanf(sampleFile, "%lf", &spectrum[i]) != EOF) // %lf tells fscanf to read a double
	{
		i++;
	}

	// Close the sample file
	fclose(sampleFile);

	// Compute the first 13 coefficients
	for(coeff = 0; coeff < 13; coeff++)
	{
		mfcc_result = GetCoefficient(spectrum, 44100, 48, 128, coeff);
		printf("%i %f\n", coeff, mfcc_result);
	}
	getchar();
	
	return 0;
}
