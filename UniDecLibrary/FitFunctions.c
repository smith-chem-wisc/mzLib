//
//  FitFunctions.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "FitFunctions.h"
#include <stdlib.h>
#include <stdio.h>
#include "ArrayIndexing.h"
#include "Convolution.h"
// Fit functions

float getfitdatspeedy(float* fitdat, const float* blur, const char* barr, const int lengthmz, const int numz, const int maxlength, const float maxint,
    const int isolength, const int* isotopepos, const float* isotopeval, const int* starttab, const int* endtab, const float* mzdist, const int speedyflag)
{
    unsigned int i, j, k;
    float* deltas = NULL;
    deltas = calloc(lengthmz, sizeof(float));
    if (isolength == 0) {
#pragma omp parallel for private (i,j), schedule(auto)
        for (i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
        {
            float temp = 0;
            for (j = 0; j < numz; j++)
            {
                temp += blur[index2D(numz, i, j)];
            }
            deltas[i] = temp;
        }
    }
    else {
        for (i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
        {
            for (j = 0; j < numz; j++)
            {
                float topval = blur[index2D(numz, i, j)];
                for (k = 0; k < isolength; k++)
                {
                    int pos = isotopepos[index3D(numz, isolength, i, j, k)];
                    float val = isotopeval[index3D(numz, isolength, i, j, k)];
                    deltas[pos] += topval * (float)val;
                }
            }
        }

    }
    if (maxlength != 0)
    {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, fitdat, speedyflag);
    }
    else
    {
        memcpy(fitdat, deltas, sizeof(float) * lengthmz);
    }

    free(deltas);
    float fitmax = 0;
    //#pragma omp parallel for private (i), schedule(auto)
    for (i = 0;i < lengthmz;i++)
    {
        if (fitdat[i] > fitmax) { fitmax = fitdat[i]; }
    }
    //#pragma omp parallel for private (i), schedule(auto)
    if (fitmax != 0)
    {
        for (i = 0; i < lengthmz; i++)
        {
            if (fitdat[i] < 0) { fitdat[i] = 0; }
            else { fitdat[i] = fitdat[i] * maxint / fitmax; }
        }
    }
    return fitmax;
}



void KillB(float* I, char* B, float intthresh, int lengthmz, int numz, 
    const int isolength, int* isotopepos, float* isotopeval)
{
    unsigned int i, j, k;
    if (isolength == 0) {
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                if (I[i] <= intthresh) { B[index2D(numz, i, j)] = 0; }
            }
        }
    }
    else
    {
        float cutoff = 0.5;
        printf("Removing species where isotope fall below %f\n", cutoff * 100);
#pragma omp parallel for private (i,j,k), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                float max = 0;
                for (k = 0; k < isolength; k++)
                {
                    float val = isotopeval[index3D(numz, isolength, i, j, k)];
                    if (val > max) { max = val; }
                }
                for (k = 0; k < isolength; k++)
                {
                    float val = isotopeval[index3D(numz, isolength, i, j, k)];
                    if (val > cutoff * max) {
                        int pos = isotopepos[index3D(numz, isolength, i, j, k)];
                        if (I[pos] <= intthresh) { B[index2D(numz, i, j)] = 0; }
                    }
                }

            }
        }
    }
}



// Fit Functions //
