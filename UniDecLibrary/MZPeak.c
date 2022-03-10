//
//  MZPeak.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "MZPeak.h"
#include <stdlib.h>
#include <stdio.h>
#include "MathUtilities.h"
#include "ArrayIndexing.h"
#include <math.h>
// MZ Peak


//Function for defining m/z peak shape. Other peak shape functions could be easily added her.
float mzpeakshape(float x, float y, float sig, int psfun)
{
    if (sig == 0) { printf("Error: mzpeakshape sigma is 0\n"); exit(103); }
    float result;
    if (psfun == 0)
    {
        result = exp(-(pow(x - y, 2)) / (2 * sig * sig));
    }
    else if (psfun == 1)
    {
        result = pow((sig / 2), 2) / (pow((x - y), 2) + pow((sig / 2), 2));
    }
    else if (psfun == 2)
    {
        if (y < x)
        {
            result = exp(-(pow(x - y, 2)) / (2 * sig * sig * 0.180337));
        }
        else
        {
            result = (sig / 2) * (sig / 2) / (pow((x - y), 2) + pow((sig / 2), 2));
        }
    }
    else
    {
        printf("Invalid Peak Function");
        exit(14);
    }
    return result;
}



void charge_scaling(float* blur, const int* nztab, const int lengthmz, const int numz)
{
    for (int i = 0; i < lengthmz; i++)
    {
        for (int j = 0; j < numz; j++)
        {
            int charge = nztab[j];
            float z = (float)charge;
            if (z != 0) { blur[index2D(numz, i, j)] /= z; }
        }
    }
    return;
}


void MakePeakShape2D(int lengthmz, int maxlength, int* starttab, int* endtab, float* dataMZ, float mzsig, int psfun, int speedyflag, float* mzdist, float* rmzdist, int makereverse)
{
    #pragma omp parallel for schedule(auto)
    for(int i=0;i<lengthmz;i++)
      {
          int start = starttab[i];
          int end = endtab[i];
          int t = 0;
          for (int j = start; j <= end; j++)
          {
              int j2 = fixk(j, lengthmz);
              mzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[i], dataMZ[j2], mzsig, psfun);
              if (makereverse == 1) { rmzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[j2], dataMZ[i], mzsig, psfun); }
          }
      }
}

void MakePeakShape1D(float* dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig, int psfun, float* mzdist, float* rmzdist, int makereverse)
{
    float binsize = dataMZ[1] - dataMZ[0];
    float newrange = threshold / binsize;
    int n;
    for (n = (int)-newrange;n < (int)newrange;n++)
    {
        mzdist[indexmod(lengthmz, 0, n)] = mzpeakshape(0, n * binsize, mzsig, psfun);
        if (makereverse == 1) { rmzdist[indexmod(lengthmz, 0, n)] = mzpeakshape(n * binsize, 0, mzsig, psfun); }
    }
    printf("\nNotice: Assuming linearized data. \n\n");
}

