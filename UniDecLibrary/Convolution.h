//
//  Convolution.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Convolution_h
#define Convolution_h

#include <stdio.h>
#include <stdlib.h>
#include "Config.h"
#include "Decon.h"
#include "fftw3.h"
#include "ArrayIndexing.h"
#include "BlurFunctions.h"
#include <math.h>
#include "Sorting.h"
#include "MathUtilities.h"

void convolve_simp(const int lengthmz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, const float* deltas, float* denom, const int speedyflag);

__declspec(dllexport) void deconvolve_baseline(const int lengthmz, const float* dataMZ, const float* dataInt, float* baseline, const float mzsig);

__declspec(dllexport) float deconvolve_iteration_speedy(const int lengthmz, const int numz, const int maxlength, const float* __restrict blur, float* __restrict blur2,
    const char* __restrict barr, const int aggressiveflag, const float* __restrict  dataInt,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval
    , const int* starttab, const int* endtab, const float* mzdist, const float* rmzdist, const int speedyflag, const int baselineflag, float* baseline,
                                  float* noise, const float mzsig, const float* dataMZ, const float filterwidth, const float psig);
__declspec(dllexport) float Reconvolve(const int lengthmz, const int numz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, const float* blur, float* newblur, const int speedyflag, const char* barr);

__declspec(dllexport) int SetStartsEnds(const Config config, const Input* inp, int* starttab, int* endtab, const float threshold);

void cconv2fast(double* a, double* b, double* c, int length);

void dd_deconv2(double* kernel_y, double* data_y, int length, double* output);

#endif /* Convolution_h */
