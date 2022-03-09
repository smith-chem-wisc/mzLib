//
//  PointSmoothing.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "PointSmoothing.h"
// Point Smoothing //


void point_smoothing(float* blur, const char* barr, const int lengthmz, const int numz, const int width)
{
    float* newblur;
    newblur = calloc(lengthmz * numz, sizeof(float));
    memcpy(newblur, blur, lengthmz * numz * sizeof(float));
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        for (int j = 0; j < numz; j++)
        {
            if (barr[index2D(numz, i, j)] == 1) {
                int low = i - width;
                if (low < 0) { low = 0; }
                int high = i + width + 1;
                if (high > lengthmz) { high = lengthmz; }

                float sum = 0;
                for (int k = low; k < high; k++)
                {
                    sum += newblur[index2D(numz, k, j)];
                }

                blur[index2D(numz, i, j)] = sum / ((float)1 + 2 * width);
            }
        }
    }
    free(newblur);
    return;
}

void point_smoothing_peak_width(const int lengthmz, const int numz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, float* blur, const int speedyflag, const char* barr)
{
    float* newblur;
    newblur = calloc(lengthmz * numz, sizeof(float));
    memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    Reconvolve(lengthmz, numz, maxlength, starttab, endtab, mzdist, newblur, blur, speedyflag, barr);
    return;
}

/*
void point_smoothing_sg(float *blur, const int lengthmz, const int numz, const int width)
{
    float *newblur;
    newblur = calloc(lengthmz*numz, sizeof(float));
    memcpy(newblur, blur, lengthmz*numz * sizeof(float));

    float *c;
    int np = width * 2 + 1;
    c = calloc(np, sizeof(float));
    savgol(c, np, width, width, 0, 2);

    //printf("Savgol: ");
    //for (int i = 0; i < np; i++){printf("%f ", c[i]);}
    //printf("\n");

    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        for (int j = 0; j < numz; j++)
        {
            int low = i - width;
            if (low < 0) { low = 0; }
            int high = i + width + 1;
            if (high > lengthmz) { high = lengthmz; }

            float sum = 0;
            for (int k = low; k < high; k++)
            {
                int cind = mod((i - k) , np);
                sum += c[cind] * newblur[index2D(numz, k, j)];
            }

            blur[index2D(numz, i, j)] = sum;
        }
    }
    free(newblur);
    return;
}
*/


// Point Smoothing //
