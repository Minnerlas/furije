#ifndef FFT_H
#define FFT_H

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define PI 3.14159265358979

void dft(double complex *l, double complex *yf, uint64_t N) {
	for(uint64_t k = 0; k<N; k++) {
		double complex t = 0;
		for(uint64_t n = 0; n<N; n++)
			t+=l[n]*cexp(-I*2*PI*k*n/N);
		yf[k] = t;
	}
}

void idft(double complex *l, double complex *yf, uint64_t N) {
	for(uint64_t k = 0; k<N; k++) {
		double complex t = 0;
		for(uint64_t n = 0; n<N; n++)
			t+=l[n]*cexp(I*2*PI*k*n/N);
		yf[k] = t/N;
	}
}

void fft(double complex *l, double complex *yf, uint64_t N) {
	//printf("%d\n", N);
	if(N%2>0) {
		printf("N must be even.\n");
		return;
	} else if(N<=4) {
		dft(l, yf, N);
		return;
	}

	double complex *lparno = malloc(sizeof(*lparno)*(N/2));
	double complex *yfparno = malloc(sizeof(*yfparno)*(N/2));
	double complex *lneparno = malloc(sizeof(*lneparno)*(N/2));
	double complex *yfneparno = malloc(sizeof(*yfneparno)*(N/2));

	for(uint64_t i = 0; i<N/2; i++)
		lparno[i] = l[2*i];
	fft(lparno, yfparno, N/2);

	for(uint64_t i = 0; i<N/2; i++)
		lneparno[i] = l[2*i+1];
	fft(lneparno, yfneparno, N/2);

	for(uint64_t i = 0; i<N; i++) {
		double complex t = cexp(-2*I*PI*i/N);
		yf[i] = yfparno[i%(N/2)] + t * yfneparno[i%(N/2)];
	}

	free(lparno);
	free(yfparno);
	free(lneparno);
	free(yfneparno);

}

void ifft(double complex *l, double complex *yf, uint64_t N) {
	//printf("%d\n", N);
	if(N%2>0) {
		printf("N must be even.\n");
		return;
	} else if(N<=4) {
		idft(l, yf, N);
		return;
	}

	double complex *lparno = malloc(sizeof(*lparno)*(N/2));
	double complex *yfparno = malloc(sizeof(*yfparno)*(N/2));
	double complex *lneparno = malloc(sizeof(*lneparno)*(N/2));
	double complex *yfneparno = malloc(sizeof(*yfneparno)*(N/2));

	for(uint64_t i = 0; i<N/2; i++)
		lparno[i] = l[2*i];
	ifft(lparno, yfparno, N/2);

	for(uint64_t i = 0; i<N/2; i++)
		lneparno[i] = l[2*i+1];
	ifft(lneparno, yfneparno, N/2);

	for(uint64_t i = 0; i<N; i++) {
		double complex t = cexp(2*I*PI*i/N);
		yf[i] = (yfparno[i%(N/2)] + t * yfneparno[i%(N/2)])/2;
	}

	free(lparno);
	free(yfparno);
	free(lneparno);
	free(yfneparno);

}

#endif /* FFT_H */
