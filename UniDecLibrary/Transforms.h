//
//  Transforms.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Transforms_h
#define Transforms_h

#include <stdio.h>
__declspec(dllexport) void InterpolateTransform(const int maaxle, const int numz, const int lengthmz, const int* nztab, float* massaxis,
                          const float adductmass, const float* dataMZ, float* massgrid, float* massaxisval, const float* blur);
__declspec(dllexport) void SmartTransform(const int maaxle, const int numz, const int lengthmz, const int* nztab, float* massaxis,
                    const float adductmass, const float* dataMZ, float* massgrid, float* massaxisval, const float* blur);
__declspec(dllexport) void IntegrateTransform(const int lengthmz, const int numz, const float* mtab, float massmax, float massmin,
                        const int maaxle, float* massaxis, float* massaxisval, const float* blur, float* massgrid); 
#endif /* Transforms_h */
