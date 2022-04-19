//
//  BlurFunctions.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "BlurFunctions.h"
#include "Config.h"
#include "MZPeak.h"
#include "Sorting.h"

//Convolution of neighborhood function with gaussian filter.
void blur_it(const int lengthmz,
    const int numz,
    const int numclose,
    const int* __restrict closeind,
    const float* __restrict closearray,
    float* __restrict newblur,
    const float* __restrict blur,
    const char* __restrict barr)
{
    if (numclose == 1)
    {
        memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    }
    else {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz * numz; i++)
        {
            float temp = 0;
            if (barr[i] == 1)
            {
                for (int k = 0; k < numclose; k++)
                {
                    if (closeind[index2D(numclose, i, k)] != -1)
                    {
                        temp += closearray[index2D(numclose, i, k)] * blur[closeind[index2D(numclose, i, k)]];
                    }
                }
            }
            newblur[i] = temp;
        }
    }
}

//Charge state smooth using a mean filter of the log
void blur_it_mean(const int lengthmz,
    const int numz,
    const int numclose,
    const int* __restrict closeind,
    float* __restrict newblur,
    const float* __restrict blur,
    const char* __restrict barr,
    const float* __restrict closearray,
    const float zerolog)
{
    if (numclose == 1)
    {
        memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    }
    else {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz * numz; i++)
        {
            float temp = 0;
            if (barr[i] == 1)
            {
                for (int k = 0; k < numclose; k++)
                {
                    float temp2 = 0;
                    if (closeind[index2D(numclose, i, k)] != -1)
                    {
                        temp2 = blur[closeind[index2D(numclose, i, k)]] * closearray[index2D(numclose, i, k)];
                    }
                    if (temp2 > 0) { temp += log(temp2); }
                    else { temp += zerolog; }
                }
                temp = exp(temp / (float)numclose);
            }
            newblur[i] = temp;
        }
    }
}

//Convolution of neighborhood function with gaussian filter.
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
    const float zerolog)
{
    int i, j, k, n;
    int numclose = zlength * mlength;
    if (numclose == 1)
    {
        memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    }
    else {
#pragma omp parallel for private (i,k,j, n), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                float temp = 0;
                if (barr[index2D(numz, i, j)] == 1)
                {
                    for (n = 0; n < mlength; n++)
                    {
                        float temp2 = 0;
                        for (k = 0; k < zlength; k++)
                        {
                            int m = index2D(mlength, k, n);
                            float temp3 = 0;
                            if (closeind[index3D(numz, numclose, i, j, m)] != -1)
                            {
                                temp3 = blur[closeind[index3D(numz, numclose, i, j, m)]] * closearray[index3D(numz, numclose, i, j, m)];
                            }
                            if (temp3 > 0) { temp2 += log(temp3); }
                            else { temp2 += zerolog; }
                        }
                        temp += exp(temp2 / (float)zlength) * mdist[n];
                    }
                }
                newblur[index2D(numz, i, j)] = temp;
            }
        }
    }
}


//Convolution of neighborhood function with gaussian filter.
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
    const float zerolog)
{
    int i, j, k, n;
    int numclose = zlength * mlength;
    if (numclose == 1)
    {
        memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    }
    else {
#pragma omp parallel for private (i,k,j, n), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                float temp = 0;
                if (barr[index2D(numz, i, j)] == 1)
                {
                    for (n = 0; n < mlength; n++)
                    {
                        float temp2 = 0;
                        for (k = 0; k < zlength; k++)
                        {
                            int m = index2D(mlength, k, n);
                            if (closeind[index3D(numz, numclose, i, j, m)] != -1)
                            {
                                temp2 += blur[closeind[index3D(numz, numclose, i, j, m)]] * zdist[k] * closearray[index3D(numz, numclose, i, j, m)];
                            }
                        }
                        if (temp2 > 0) { temp += log(temp2); }// / (float)mlength);}
                        else { temp += zerolog; }
                    }
                    temp = exp(temp / (float)mlength);
                }
                newblur[index2D(numz, i, j)] = temp;
            }
        }
    }
}

void midblur_baseline(float* baseline, const int lengthmz, const float* dataMZ, const float mzsig, int mult)
{
    if (mult == 0)
    {
        mult = lengthmz / 400;
    }
    int window = 25;

    float* temp = NULL;
    temp = calloc(lengthmz, sizeof(float));
    memcpy(temp, baseline, sizeof(float) * lengthmz);
    int i, j;
#pragma omp parallel for private(i,j), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        float val = 0;
        float len = window * 2;
        float* med = NULL;
        med = calloc(len, sizeof(float));
        int index = 0;
        for (j = -window; j < window; j++)
        {
            int k = i + j * mult;
            float newval = 0;
            if (k >= 0 && k < lengthmz) {
                newval = temp[k];
            }
            if (k < 0)
            {
                newval = temp[-k];
            }
            if (k >= lengthmz)
            {
                newval = temp[2 * lengthmz - k];
            }
            med[index] = newval;
            index++;
        }
        qsort(med, len, sizeof(float), compare_function);
        index = 0;
        for (j = 0; j < window; j++)
        {
            val += med[j];
            index++;
        }
        if (index != 0) {
            baseline[i] = val / ((float)index);
        }
        else { baseline[i] = 0; }
        free(med);
    }
    free(temp);
}

void blur_noise(float* noise, int lengthmz)
{
    float* temp = NULL;
    temp = calloc(lengthmz, sizeof(float));
    memcpy(temp, noise, sizeof(float) * lengthmz);
    int i, j;
    float filter[5] = { -0.1,-0.4,1,-0.4,-0.1 };
#pragma omp parallel for private (i,j), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        float val = 0;
        for (j = -2; j <= 2; j++)
        {
            int k = i + j;
            float newval = 0;
            if (k >= 0 && k < lengthmz) {
                newval = temp[k];
            }
            if (k < 0)
            {
                newval = temp[-k];
            }
            if (k >= lengthmz)
            {
                newval = temp[2 * lengthmz - k];
            }
            val += newval * filter[j + 2];
        }
        noise[i] = val;
    }
    free(temp);
}


// Blur functions

void blur_baseline(float* baseline, const int lengthmz, const float* dataMZ, const float mzsig, int mult, int filterwidth)
{
    int mulin = mult;
    float* temp = NULL;
    temp = calloc(lengthmz, sizeof(float));
    memcpy(temp, baseline, sizeof(float) * lengthmz);
    int i, j;
#pragma omp parallel for private (i,j), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        float mzdiff = 0;
        if (i > 0) {
            mzdiff = dataMZ[i] - dataMZ[i - 1];
        }
        else {
            mzdiff = dataMZ[i + 1] - dataMZ[i];
        }
        if (mulin == 0 && mzdiff > 0)
        {
            //mult = lengthmz / 500;
            mult = (int)(2 * mzsig / mzdiff);
        }
        if (mult < 1) { mult = 1; }

        float val = 0;
        int window = filterwidth;
        for (j = -window; j < window; j++)
        {
            int k = i + j * (mult);
            float newval = 0;
            if (k >= 0 && k < lengthmz) {
                newval = temp[k];
            }
            if (k < 0)
            {
                newval = temp[-k];
            }
            if (k >= lengthmz)
            {
                newval = temp[2 * lengthmz - k];
            }
            //if(newval>0){
            val += newval;
            //}
        }
        baseline[i] = val / ((float)window * 2 + 1);
    }
    free(temp);
}


void MakeSparseBlur(const int numclose, char* barr, const int* closezind,
    const int* closemind, const float* mtab, const int* nztab,
    const float* dataMZ, int* closeind, float* closeval, float* closearray, const Config config)
{
    int lengthmz = config.lengthmz;
    int numz = config.numz;
    float molig = config.molig;
    float adductmass = config.adductmass;

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        for (int j = 0; j < numz; j++)
        {
            if (barr[index2D(numz, i, j)] == 1)
            {
                int num = 0;

                //Reset the threshold if it is zero
                float mzsig = config.mzsig;
                if (mzsig == 0) {
                    int i1 = i - 1;
                    int i2 = i + 1;
                    if (i >= lengthmz - 1) { i2 = i; }
                    if (i == 0) { i1 = i; }
                    mzsig = 2 * fabs(dataMZ[i2] - dataMZ[i1]);
                    if (mzsig > config.massbins || mzsig == 0) { mzsig = config.massbins * 2; }
                }
                float newthreshold = mzsig * 2;

                for (int k = 0; k < numclose; k++)
                {
                    //Find the z index and test if it is within the charge range
                    int indz = (int)(j + closezind[k]);
                    if (indz < 0 || indz >= numz || (nztab[j] + closezind[k]) == 0)
                    {
                        closeind[index3D(numz, numclose, i, j, k)] = -1;
                        closearray[index3D(numz, numclose, i, j, k)] = 0;
                    }
                    else
                    {
                        //Find the nearest m/z value and test if it is close enough to the predicted one and within appropriate ranges
                        float point = (float)((mtab[index2D(numz, i, j)] + closemind[k] * molig + adductmass * (float)(nztab[j] + closezind[k])) / (float)(nztab[j] + closezind[k]));
                        if (point<dataMZ[0] - newthreshold || point>dataMZ[lengthmz - 1] + newthreshold)
                        {
                            closeind[index3D(numz, numclose, i, j, k)] = -1;
                            closearray[index3D(numz, numclose, i, j, k)] = 0;
                        }
                        else {
                            int ind = nearfast(dataMZ, point, lengthmz);
                            float closepoint = dataMZ[ind];
                            int newind = index2D(numz, ind, indz);
                            if (barr[newind] == 1 && fabs(point - closepoint) < newthreshold) {
                                closeind[index3D(numz, numclose, i, j, k)] = newind;
                                closearray[index3D(numz, numclose, i, j, k)] = closeval[k] * mzpeakshape(point, closepoint, mzsig, config.psfun);
                                num += 1;
                            }
                            else {
                                closeind[index3D(numz, numclose, i, j, k)] = -1;
                                closearray[index3D(numz, numclose, i, j, k)] = 0;
                            }
                            //printf("%d %d %d %f %f %d\n", i, j, k, point, closepoint, closeind[index3D(numz, numclose, i, j, k)]);
                        }

                    }
                }
                if (num < 2 && config.isotopemode == 0) { barr[index2D(numz, i, j)] = 0; }// printf("%d %d \n", i, j);}
            }
            else
            {
                for (int k = 0; k < numclose; k++)
                {
                    closeind[index3D(numz, numclose, i, j, k)] = -1;
                    closearray[index3D(numz, numclose, i, j, k)] = 0;
                }
            }
        }
    }
}


