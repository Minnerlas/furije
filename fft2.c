#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "fft.h"

int main() {
	int N = 0, fsamp = 0, cplx;
	FILE *ul = fopen("izlaz.dat", "r");

	fscanf(ul, "%d %d %d\n", &N, &fsamp, &cplx);
	printf("%d %d %d\n", N, fsamp, cplx);

	double complex *l = malloc(sizeof(*l)*N);
	double complex *yf = malloc(sizeof(*yf)*N);

	for(int i = 0; i<N; i++) {
		double r, im;
		fscanf(ul, "%lf %lf\n", &r, &im);
		l[i] = r + im * I;
	}

	fft2(l, yf, N);

	// for(int i = 0; i<10; i++)
	// 	printf("%lf %lf\n", creal(l[i]), cimag(l[i]));

	free(l);
	fclose(ul);

	FILE *iz = fopen("fft2.dat", "w");

	fprintf(iz, "%d %d %d", N, fsamp, 1);
	for(int i = 0; i<N; i++)
		fprintf(iz, "\n%.14lf %.14lf", creal(yf[i]), cimag(yf[i]));

	fclose(iz);
	free(yf);
}
