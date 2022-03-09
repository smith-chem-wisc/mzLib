//
//  ArrayIndexing.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef ArrayIndexing_h
#define ArrayIndexing_h

#include <stdio.h>
inline int index2D(const int ncols, const int r, const int c);
inline int indexmod(int length, int r, int c);
inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d);
float max1d(float* blur, int lengthmz);
inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d);
float max1d(float* blur, int lengthmz);
float min1d(float* blur, int lengthmz);
void ApplyCutoff1D(float* array, float cutoff, int lengthmz);

#endif /* ArrayIndexing_h */
