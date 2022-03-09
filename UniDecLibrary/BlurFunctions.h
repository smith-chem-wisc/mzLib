//
//  BlurFunctions.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef BlurFunctions_h
#define BlurFunctions_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void blur_it(const int lengthmz,
    const int numz,
    const int numclose,
    const int* __restrict closeind,
    const float* __restrict closearray,
    float* __restrict newblur,
    const float* __restrict blur,
             const char* __restrict barr); 

void blur_it_mean(const int lengthmz,
    const int numz,
    const int numclose,
    const int* __restrict closeind,
    float* __restrict newblur,
    const float* __restrict blur,
    const char* __restrict barr,
    const float* __restrict closearray,
                  const float zerolog);

void blur_it_hybrid1(const int lengthmz,
    const int numz,
    const int zlength,
    const int mlength,
    const int* __restrict closeind,
    const int* __restrict closemind,
    const int* __restrict closezind,
    const float* __restrict mdist,
    const float* __restrict zdist,
    float* __restrict newblur,
    const float* __restrict blur,
    const char* __restrict barr,
    const float* __restrict closearray,
                     const float zerolog);

void blur_it_hybrid2(const int lengthmz,
    const int numz,
    const int zlength,
    const int mlength,
    const int* __restrict closeind,
    const int* __restrict closemind,
    const int* __restrict closezind,
    const float* __restrict mdist,
    const float* __restrict zdist,
    float* __restrict newblur,
    const float* __restrict blur,
    const char* __restrict barr,
    const float* __restrict closearray,
                     const float zerolog);


void midblur_baseline(float* baseline, const int lengthmz, const float* dataMZ, const float mzsig, int mult);
void blur_noise(float* noise, int lengthmz);
void blur_baseline(float* baseline, const int lengthmz, const float* dataMZ, const float mzsig, int mult, int filterwidth);
void MakeSparseBlur(const int numclose, char* barr, const int* closezind,
    const int* closemind, const float* mtab, const int* nztab,
                    const float* dataMZ, int* closeind, float* closeval, float* closearray, const Config config);
void midblur_baseline(float* baseline, const int lengthmz, const float* dataMZ, const float mzsig, int mult);

#endif /* BlurFunctions_h */
