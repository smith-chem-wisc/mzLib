//
//  ArrayIndexing.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "ArrayIndexing.h"
// Array indexing functions:
inline int index2D(const int ncols, const int r, const int c)
{
    return r * ncols + c;
}

inline int indexmod(int length, int r, int c)
{
    return mod((c - r), length);
}

inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d)
{
    return r * ncols * nrows + c * nrows + d;
}


float max1d(float* blur, int lengthmz) {
    float blurmax = blur[0];
    unsigned int i;
    for (i = 0; i < lengthmz; i++)
    {
        if (blur[i] > blurmax)
        {
            blurmax = blur[i];
        }
    }
    return blurmax;
}

float min1d(float* blur, int lengthmz) {
    float blurmin = blur[0];
    unsigned int i;
    for (i = 0; i < lengthmz; i++)
    {
        if (blur[i] < blurmin)
        {
            blurmin = blur[i];
        }
    }
    return blurmin;
}


void ApplyCutoff1D(float* array, float cutoff, int lengthmz)
{
    unsigned int i;
    //#pragma omp parallel for private (i,j), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        if (array[i] < cutoff) { array[i] = 0; }
    }
}

// Array Indexing

