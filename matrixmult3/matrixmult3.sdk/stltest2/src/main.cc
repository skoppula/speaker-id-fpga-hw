/*
 * Empty C++ Application
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main()
{
	return 0;
}
