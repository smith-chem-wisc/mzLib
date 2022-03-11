//
//  MathUtilities.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef MathUtilities_h
#define MathUtilities_h

#include <stdio.h>
#include "fftw3.h"
#include "ArrayIndexing.h"
#include <math.h>

float secderndis(float m, float s, float x);
int compare_function(const void* a, const void* b);

void sum_deltas(const int lengthmz, const int numz, const float* __restrict blur, const char* __restrict barr,
                const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, float* deltas);

void apply_ratios(const int lengthmz, const int numz, const float* __restrict blur, const char* __restrict barr,
                  const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, const float* __restrict denom, float* blur2);

void complex_mult(fftw_complex* A, fftw_complex* B, fftw_complex* product_ft, int complen);

__declspec(dllexport) float Average(const int length, const float* xarray);
int mod(int a, int b); 
float Max(const float* blur, const int length);
float Sum(const float* blur, int length);
int twopow(int num);
float nativecharge(float mass, float fudge);
float StdDev(int length, float* xarray, float wmean);
float ndis(float x, float y, float sig);
float clip(float x, float cutoff);
void ignorezeros(char* barr, const float* dataInt, const int lengthmz, const int numz);
inline int fixk(int k, int lengthmz)
{
    k = abs(k);
    if (k >= lengthmz) { k = 2 * lengthmz - k - 2; }
    //if (k < 0) { k = 0; }
    return k;
}
#endif /* MathUtilities_h */
