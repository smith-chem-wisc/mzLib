//
//  PointSmoothing.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef PointSmoothing_h
#define PointSmoothing_h

#include <stdio.h>
void point_smoothing(float* blur, const char* barr, const int lengthmz, const int numz, const int width);
void point_smoothing_peak_width(const int lengthmz, const int numz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, float* blur, const int speedyflag, const char* barr);


#endif /* PointSmoothing_h */
