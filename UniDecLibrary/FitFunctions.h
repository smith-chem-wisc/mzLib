//
//  FitFunctions.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef FitFunctions_h
#define FitFunctions_h

#include <stdio.h>
float getfitdatspeedy(float* fitdat, const float* blur, const char* barr, const int lengthmz, const int numz, const int maxlength, const float maxint,
                      const int isolength, const int* isotopepos, const float* isotopeval, const int* starttab, const int* endtab, const float* mzdist, const int speedyflag);
void KillB(float* I, char* B, float intthresh, int lengthmz, int numz, const int isolength, int* isotopepos, float* isotopeval);

#endif /* FitFunctions_h */
