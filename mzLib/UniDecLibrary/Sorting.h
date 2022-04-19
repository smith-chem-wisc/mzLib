//
//  Sorting.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Sorting_h
#define Sorting_h

#include <stdio.h>
__declspec(dllexport) int nearunsorted(float* testmasses, float point, int lengthtest);
__declspec(dllexport) int neartest(float* testmasses, float point, int lengthtest, float cutoff);
__declspec(dllexport) int nearfast(const float* dataMZ, const float point, const int numdat);
__declspec(dllexport) int nearfast_test(const float* dataMZ, const float point, const int numdat, float cutoff);
__declspec(dllexport) int nearfast_d(const double* dataMZ, const double point, const int numdat);
#endif /* Sorting_h */
