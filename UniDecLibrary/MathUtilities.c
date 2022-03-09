//
//  MathUtilities.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "MathUtilities.h"

// General math //

//Second Derivative of a Gaussian
float secderndis(float m, float s, float x)
{
    if (s == 0) { return 0; }
    return s / 2 * (2 * exp(-pow(m - x, 2) / s)) / s - (4 * exp(-pow(m - x, 2) / s) * pow(m - x, 2)) / pow(s, 2);
}


int compare_function(const void* a, const void* b)
{
    float* x = (float*)a;
    float* y = (float*)b;
    if (*x < *y) { return -1; }
    else if (*x > * y) { return 1; }
    return 0;
}


inline int fixk(int k, int lengthmz)
{
    k = abs(k);
    if (k >= lengthmz) { k = 2 * lengthmz - k - 2; }
    //if (k < 0) { k = 0; }
    return k;
}


void sum_deltas(const int lengthmz, const int numz, const float* __restrict blur, const char* __restrict barr,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, float* deltas)
{
    int i, j, k;
    float temp;
    //Collapse the grid into a 1D array of delta function values
    if (isolength == 0) {
#pragma omp parallel for private (i,j), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            temp = 0;
            for (j = 0; j < numz; j++)
            {
                if (barr[index2D(numz, i, j)] == 1)
                {
                    temp += blur[index2D(numz, i, j)];
                }
            }
            deltas[i] = temp;
        }
    }
    else {
        //#pragma omp parallel for private (i,j,k), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                if (barr[index2D(numz, i, j)] == 1) {
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
    }
}

void apply_ratios(const int lengthmz, const int numz, const float* __restrict blur, const char* __restrict barr,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, const float* __restrict denom, float* blur2)
{
    int i, j, k;
    if (isolength == 0)
    {
#pragma omp parallel for private (i,j), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                if (barr[index2D(numz, i, j)] == 1)
                {
                    blur2[index2D(numz, i, j)] = denom[i] * blur[index2D(numz, i, j)];
                }
                else
                {
                    blur2[index2D(numz, i, j)] = 0;
                }
            }
        }
    }
    else
    {
#pragma omp parallel for private (i,j,k), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                if (barr[index2D(numz, i, j)] == 1)
                {
                    float temp = 0;
                    for (k = 0; k < isolength; k++)
                    {
                        int pos = isotopepos[index3D(numz, isolength, i, j, k)];
                        float val = isotopeval[index3D(numz, isolength, i, j, k)];
                        temp += val * (denom[pos]);
                    }
                    blur2[index2D(numz, i, j)] = (temp)*blur[index2D(numz, i, j)];
                }
                else
                {
                    blur2[index2D(numz, i, j)] = 0;
                }
            }
        }
    }
}


void complex_mult(fftw_complex* A, fftw_complex* B, fftw_complex* product_ft, int complen) {
    // A * B = (ac - bd) + i(ad + bc)
    for (int j = 0; j < complen; j++) {
        product_ft[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
        product_ft[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
    }
}


//Calculate Average
float Average(const int length, const float* xarray)
{
    float temp1 = 0;
    float temp2 = (float)length;
    for (int i = 0; i < length; i++)
    {
        temp1 += xarray[i];
    }
    if (temp2 == 0) { return 0; }
    return temp1 / temp2;
}

//Actual Modulus operator rather than Remainder operator %
int mod(int a, int b) { int r = a % b; return r < 0 ? r + b : r; }

float Max(const float* blur, const int length) {
    float blurmax = 0;
    for (int i = 0; i < length; i++)
    {
        if (blur[i] > blurmax)
        {
            blurmax = blur[i];
        }
    }
    return blurmax;
}

float Sum(const float* blur, int length) {
    float sum = 0;
    for (int i = 0; i < length; i++)
    {
        sum += blur[i];
    }
    return sum;
}


//Finds nearest factor of two for optimizing FFTs
int twopow(int num) {
    int n, val;
    n = 0;
    val = 1;
    while (n < 100 && val < num) {
        n++;
        val = pow(2, n);
    }
    return val;
}

//Average native charge state from Champ
float nativecharge(float mass, float fudge)
{
    return 0.0467 * pow(mass, 0.533) + fudge;
}


//Calculate Standard Deviation
float StdDev(int length, float* xarray, float wmean)
{
    float temp1 = 0;
    float temp2 = 0;
    int i;
    for (i = 0;i < length;i++)
    {
        temp1 += pow(xarray[i] - wmean, 2);
        temp2 += 1;
    }
    if (temp2 == 0) { return 0; }
    return sqrt(temp1 / temp2);
}

float ndis(float x, float y, float sig)
{
    if (sig == 0) { return 0; }
    return 1 / (sig * 2.50663) * exp(-(pow(x - y, 2.)) / (2. * sig * sig));
}

float clip(float x, float cutoff)
{
    if (x > cutoff) { return x; }
    else { return 0; }
}


void ignorezeros(char* barr, const float* dataInt, const int lengthmz, const int numz)
{
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++)
    {
        float val = dataInt[i];
        if (val == 0) {
            for (int j = 0; j < numz; j++)
            {
                barr[index2D(numz, i, j)] = 0;
            }
        }
    }
}
