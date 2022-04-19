//
//  MZPeak.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef MZPeak_h
#define MZPeak_h

#include <stdio.h>
float mzpeakshape(float x, float y, float sig, int psfun);
void charge_scaling(float* blur, const int* nztab, const int lengthmz, const int numz);
__declspec(dllexport) void MakePeakShape2D(int lengthmz, int maxlength, int* starttab, int* endtab, float* dataMZ, float mzsig, int psfun, int speedyflag, float* mzdist, float* rmzdist, int makereverse);
__declspec(dllexport) void MakePeakShape1D(float* dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig, int psfun, float* mzdist, float* rmzdist, int makereverse);

#endif /* MZPeak_h */
