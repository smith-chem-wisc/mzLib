//
//  Decon.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Decon_h
#define Decon_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Config.h"
#include "Input.h"

__declspec(dllexport) typedef struct Decon Decon;

struct Decon {
    float* fitdat;
    float* baseline;
    float* noise;
    float* massgrid;
    float* massaxis;
    float* massaxisval;
    float* blur;
    float* newblur;
    float* peakx;
    float* peaky;
    float* dscores;
    float error;
    float rsquared;
    int iterations;
    float uniscore;
    float conv;
    float threshold;
    int mlen;
    int plen;
};

__declspec(dllexport) Decon SetupDecon(void);

__declspec(dllexport) void FreeDecon(Decon decon);

#endif /* Decon_h */
