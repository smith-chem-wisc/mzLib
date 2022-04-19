#include "TestFFTW3.h"

int ExecuteTestFFT(void) {
	int N = 100; 
	fftw_complex* in, * out; 
	fftw_plan fft_plan; 
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); 
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); 
	fft_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 

	fftw_execute(fft_plan); 

	fftw_destroy_plan(fft_plan); 
	fftw_free(in); 
	fftw_free(out); 
	
	return 0; 

}