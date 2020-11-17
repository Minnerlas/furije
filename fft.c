#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#define PI 3.14159265358979

void dft(double complex *l, double complex *yf, int N) {
	for(int k = 0; k<N; k++) {
		double complex t = 0;
		for(int n = 0; n<N; n++)
			t+=l[n]*cexp(-I*2*PI*k*n/N);
		yf[k] = t;
	}
}

void fft(double complex *l, double complex *yf, int N) {
	//printf("%d\n", N);
	if(N<=4) {
		dft(l, yf, N);
		return;
	}

	double complex *lparno = malloc(sizeof(*lparno)*(N/2));
	double complex *yfparno = malloc(sizeof(*yfparno)*(N/2));
	double complex *lneparno = malloc(sizeof(*lneparno)*(N/2));
	double complex *yfneparno = malloc(sizeof(*yfneparno)*(N/2));

	for(int i = 0; i<N/2; i++)
		lparno[i] = l[2*i];
	fft(lparno, yfparno, N/2);

	for(int i = 0; i<N/2; i++)
		lneparno[i] = l[2*i+1];
	fft(lneparno, yfneparno, N/2);

	for(int i = 0; i<N; i++) {
		double complex t = cexp(-2*I*PI*i/N);
		yf[i] = yfparno[i%(N/2)] + t * yfneparno[i%(N/2)];
	}

	free(lparno);
	free(yfparno);
	free(lneparno);
	free(yfneparno);

}

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

	fft(l, yf, N);

	// for(int i = 0; i<10; i++)
	// 	printf("%lf %lf\n", creal(l[i]), cimag(l[i]));

	free(l);
	fclose(ul);

	FILE *iz = fopen("fft_c.dat", "w");

	fprintf(iz, "%d %d %d", N, fsamp, 1);
	for(int i = 0; i<N; i++)
		fprintf(iz, "\n%.14lf %.14lf", creal(yf[i]), cimag(yf[i]));

	fclose(iz);
	free(yf);
}
