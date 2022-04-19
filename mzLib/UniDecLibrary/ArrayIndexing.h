//
//  ArrayIndexing.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef ArrayIndexing_h
#define ArrayIndexing_h

#include <stdio.h>
#include <stdlib.h>
#include "MathUtilities.h"
#include "ArrayIndexing.h"

__declspec(dllexport) float max1d(float* blur, int lengthmz);

float min1d(float* blur, int lengthmz);

void ApplyCutoff1D(float* array, float cutoff, int lengthmz);

inline int index2D(const int ncols, const int r, const int c) {
    return r * ncols + c;
}; 
inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d)
{
    return r * ncols * nrows + c * nrows + d;
}
inline int indexmod(int length, int r, int c)
{
    int a = c - r; 
    int b = length; 
    int result = a % b; 
    return result < 0 ? result + b : result;
};

#endif /* ArrayIndexing_h */
