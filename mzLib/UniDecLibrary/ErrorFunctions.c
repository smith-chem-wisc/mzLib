//
//  ErrorFunctions.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "ErrorFunctions.h"
#include "MZPeak.h"
#include "ArrayIndexing.h"
#include "MathUtilities.h"
#include "ErrorFunctions.h"


// Error calculations

float errfunspeedy(Config config, Decon decon, const char* barr,
    const float* dataInt, const int maxlength, const int* isotopepos, 
    const float* isotopeval, const int* starttab, const int* endtab, 
    const float* mzdist, float rsquared, float* fitdat, float* blur, float* baseline)
{
    //Get max intensity
    float maxint = 0;
    for (int i = 0; i < config.lengthmz; i++)
    {
        if (dataInt[i] > maxint) { maxint = dataInt[i]; }
    }

    float fitmax = getfitdatspeedy(fitdat, blur, barr, config.lengthmz, config.numz, maxlength,
        maxint, config.isolength, isotopepos, isotopeval, starttab, endtab, mzdist, config.speedyflag);

    if (config.baselineflag == 1)
    {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < config.lengthmz; i++) {
            fitdat[i] += baseline[i];// +decon.noise[i];
            //decon.fitdat[i] = decon.noise[i]+0.1;
        }
    }
    ApplyCutoff1D(fitdat, 0, config.lengthmz);

    float fitmean = Average(config.lengthmz, dataInt);

    float error = 0;
    float sstot = 0;
    for (int i = 0;i < config.lengthmz;i++)
    {
        error += pow((fitdat[i] - dataInt[i]), 2);
        sstot += pow((dataInt[i] - fitmean), 2);
    }

    //Calculate R-squared
    if (sstot != 0) { decon.rsquared = 1 - (error / sstot); }

    return error;
}

