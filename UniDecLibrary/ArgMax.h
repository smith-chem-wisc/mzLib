//
//  ArgMax.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef ArgMax_h
#define ArgMax_h

#include <stdio.h>

int argmax(float* blur, int lengthmz); 
void softargmax_transposed(float* blur, const int lengthmz, const int numz, const float beta, const char* barr, const int maxlength,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, const int speedyflag
                           , int* starttab, int* endtab, float* mzdist, const float mzsig);
void softargmax_everything(float* blur, const int lengthmz, const int numz, const float beta);
void softargmax(float* blur, const int lengthmz, const int numz, const float beta);

#endif /* ArgMax_h */
