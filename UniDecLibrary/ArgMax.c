//
//  ArgMax.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "ArgMax.h"


// Softargmax //

int argmax(float* blur, int lengthmz)
{
    float max = blur[0];
    int pos = 0;
    unsigned int i;
    for (i = 0; i < lengthmz; i++)
    {
        if (blur[i] > max)
        {
            max = blur[i];
            pos = i;
        }
    }
    return pos;
}


void softargmax_transposed(float* blur, const int lengthmz, const int numz, const float beta, const char* barr, const int maxlength,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, const int speedyflag
    , int* starttab, int* endtab, float* mzdist, const float mzsig)
{
    float* newblur, * deltas, * deltas2, * denom, * denom2;
    newblur = calloc(lengthmz * numz, sizeof(float));
    deltas = calloc(lengthmz, sizeof(float));
    deltas2 = calloc(lengthmz, sizeof(float));
    denom = calloc(lengthmz, sizeof(float));
    denom2 = calloc(lengthmz, sizeof(float));
    memcpy(newblur, blur, lengthmz * numz * sizeof(float));

    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas);

    if (mzsig != 0)//Convolve with peak shape
    {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
    }
    else { memcpy(denom, deltas, sizeof(float) * lengthmz); }

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz * numz; i++)
    {
        blur[i] = exp(beta * newblur[i]) - 1;
    }

    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas2);

    if (mzsig != 0)//Convolve with peak shape
    {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas2, denom2, speedyflag);
    }
    else { memcpy(denom2, deltas, sizeof(float) * lengthmz); }

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        float factor = 0;
        if (denom2[i] != 0) { factor = denom[i] / denom2[i]; };
        for (int j = 0; j < numz; j++)
        {
            //blur[index2D(numz, i, j)] -= 1;
            blur[index2D(numz, i, j)] *= factor;
        }
        //printf("%f %f\n", sum1, sum3);
    }
    free(newblur);
    free(deltas);
    free(deltas2);
    free(denom);
    free(denom2);
    return;
}

void softargmax_everything(float* blur, const int lengthmz, const int numz, const float beta)
{
    float* newblur;
    newblur = calloc(lengthmz * numz, sizeof(float));
    memcpy(newblur, blur, lengthmz * numz * sizeof(float));
    //float max1 = 0;
    //float max2 = 0;
    float sum2 = 0;
    float sum1 = 0;
    float min2 = 100000000000000000;
    //#pragma omp parallel for schedule(auto), reduction(min:min2), reduction(+:sum1), reduction(+:sum2)
    for (int i = 0; i < lengthmz * numz; i++)
    {
        float d = newblur[i];
        //if (d > max1) { max1 = d; }
        float e = exp(beta * d);
        //float e = pow(d, -beta);
        //if (e > max2) { max2 = e; }
        blur[i] = e;
        if (e < min2) { min2 = e; }
        sum1 += d;
        sum2 += e;
    }
    float factor = 0;
    //float denom = (max2 - min2);
    //if (denom != 0) { factor = max1 / denom; };
    float denom = (sum2 - min2 * numz * lengthmz);
    if (denom != 0) { factor = sum1 / denom; };
    //if (sum2 != 0) { factor = sum1 / sum2; };
    if (factor > 0) {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz * numz; i++)
        {
            blur[i] -= min2;
            blur[i] *= factor;
        }
    }
    free(newblur);
    return;
}

void softargmax(float* blur, const int lengthmz, const int numz, const float beta)
{
    if (beta < 0) {
        //softargmax_transposed(blur, lengthmz, numz, fabs(beta));
        softargmax_everything(blur, lengthmz, numz, fabs(beta));
        return;
    }

    float* newblur;
    newblur = calloc(lengthmz * numz, sizeof(float));
    memcpy(newblur, blur, lengthmz * numz * sizeof(float));
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        float sum2 = 0;
        float sum1 = 0;
        float factor = 0;
        float min2 = 1000000000000;

        for (int j = 0; j < numz; j++)
        {
            float d = newblur[index2D(numz, i, j)];
            sum1 += d;

            float e = exp(beta * d);

            if (e < min2) { min2 = e; }

            blur[index2D(numz, i, j)] = e;
            sum2 += e;
        }

        float denom = (sum2 - min2 * numz);
        if (denom != 0) { factor = sum1 / denom; };

        if (factor > 0) {
            for (int j = 0; j < numz; j++)
            {
                blur[index2D(numz, i, j)] -= min2;
                blur[index2D(numz, i, j)] *= factor;
            }
        }
        else {
            for (int j = 0; j < numz; j++)
            {
                blur[index2D(numz, i, j)] = 0;
            }
        }
    }
    free(newblur);
    return;
}



// Softargmax //

