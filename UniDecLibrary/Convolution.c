//
//  Convolution.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Convolution.h"

// Convolution

void convolve_simp(const int lengthmz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, const float* deltas, float* denom, const int speedyflag)
{
    if (speedyflag == 0) {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++)
        {
            float cv = 0;
            for (int k = starttab[i]; k <= endtab[i]; k++)
            {
                int k2 = fixk(k, lengthmz);
                int start = starttab[k2];
                cv += deltas[k2] * mzdist[index2D(maxlength, k2, i - start)];
            }
            denom[i] = cv;
        }
    }
    else {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++)
        {
            float cv = 0;
            for (int k = starttab[i]; k <= endtab[i]; k++)
            {
                cv += deltas[k] * mzdist[indexmod(lengthmz, k, i)];
            }
            denom[i] = cv;
        }
    }
}


void deconvolve_baseline(const int lengthmz, const float* dataMZ, const float* dataInt, float* baseline, const float mzsig)
{
    float* denom = NULL;
    denom = calloc(lengthmz, sizeof(float));

    midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 0);
    midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 5);

    memcpy(denom, baseline, sizeof(float) * lengthmz);
    int i;
#pragma omp parallel for private(i), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
    }

    midblur_baseline(denom, lengthmz, dataMZ, mzsig, 0);
    midblur_baseline(denom, lengthmz, dataMZ, mzsig, 5);
#pragma omp parallel for private(i), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        baseline[i] = baseline[i] * (denom[i]);
    }
    free(denom);
}


float deconvolve_iteration_speedy(const int lengthmz, const int numz, const int maxlength, const float* __restrict blur, float* __restrict blur2,
    const char* __restrict barr, const int aggressiveflag, const float* __restrict  dataInt,
    const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval
    , const int* starttab, const int* endtab, const float* mzdist, const float* rmzdist, const int speedyflag, const int baselineflag, float* baseline,
    float* noise, const float mzsig, const float* dataMZ, const float filterwidth, const float psig)
{
    unsigned int i, j, k;
    float* deltas = NULL, * denom = NULL;
    deltas = calloc(lengthmz, sizeof(float));
    denom = calloc(lengthmz, sizeof(float));

    if (aggressiveflag == 1 && mzsig != 0) {
        blur_baseline(baseline, lengthmz, dataMZ, fabs(mzsig), 0, filterwidth);
        //blur_baseline(baseline, lengthmz, 10);
        //blur_noise(noise, lengthmz);
    }
    //printf("1\n");
    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas);
    //printf("2\n");
    if (mzsig != 0 && psig >= 0)
    {
        //Convolve with peak shape
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
    }
    else
    {
        memcpy(denom, deltas, sizeof(float) * lengthmz);
    }
    //printf("3\n");
    if (aggressiveflag == 1)
    {
#pragma omp parallel for private(i), schedule(auto)
        for (i = 0;i < lengthmz;i++) {
            denom[i] += baseline[i];// +noise[i]);
        }
    }
    //printf("4\n");
    //Calculate Ratio
#pragma omp parallel for private(i), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
    }
    //printf("5\n");
    if (mzsig < 0)
    {
        //Real Richardson-Lucy Second Convolution
        convolve_simp(lengthmz, maxlength, starttab, endtab, rmzdist, denom, deltas, speedyflag);
        memcpy(denom, deltas, sizeof(float) * lengthmz);
    }

    //Multiply Ratio by prior
    apply_ratios(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, denom, blur2);

    if (aggressiveflag == 1)
    {
        blur_baseline(denom, lengthmz, dataMZ, fabs(mzsig), 0, filterwidth);

#pragma omp parallel for private(i), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            baseline[i] = baseline[i] * (denom[i]);
        }
    }

    free(deltas);
    free(denom);
    return 0;
}


float Reconvolve(const int lengthmz, const int numz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, const float* blur, float* newblur, const int speedyflag, const char* barr)
{
    float newblurmax = 0;
    if (speedyflag == 0) {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++)
        {
            for (int j = 0; j < numz; j++)
            {

                float cv = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int k = starttab[i]; k <= endtab[i]; k++)
                    {
                        int k2 = fixk(k, lengthmz);
                        if (blur[index2D(numz, k2, j)] != 0)
                        {
                            int start = starttab[k2];
                            cv += blur[index2D(numz, k2, j)] * mzdist[index2D(maxlength, k2, i - start)];
                        }
                    }
                }
                newblur[index2D(numz, i, j)] = cv;
                if (cv > newblurmax)
                {
                    newblurmax = cv;
                }
            }
        }
    }
    else {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++)
        {
            for (int j = 0; j < numz; j++)
            {

                float cv = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int k = starttab[i]; k <= endtab[i]; k++)
                    {
                        if (blur[index2D(numz, k, j)] != 0)
                        {
                            cv += blur[index2D(numz, k, j)] * mzdist[indexmod(lengthmz, k, i)];
                        }
                    }
                }
                newblur[index2D(numz, i, j)] = cv;
                if (cv > newblurmax)
                {
                    newblurmax = cv;
                }
            }
        }

    }
    return newblurmax;
}

//Sets the maxlength parameter and the start and end values for the m/z peak shape convolution
//Convolution uses a reflection for the edges, so some care needs to be taken when things are over the edge.
int SetStartsEnds(const Config config, const Input* inp, int* starttab, int* endtab, const float threshold) {
    int maxlength = 1;
    for (int i = 0; i < config.lengthmz; i++)
    {
        float point = inp->dataMZ[i] - threshold;
        int start, end;
        if (point < inp->dataMZ[0] && config.speedyflag == 0) {
            //start = (int)((point - inp->dataMZ[0]) / (inp->dataMZ[1] - inp->dataMZ[0]));
            start = 0 - nearfast(inp->dataMZ, (float)2 * inp->dataMZ[0] - point, config.lengthmz);
        }
        else {
            start = nearfast(inp->dataMZ, point, config.lengthmz);
        }
        starttab[i] = start;

        point = inp->dataMZ[i] + threshold;
        if (point > inp->dataMZ[config.lengthmz - 1] && config.speedyflag == 0) {
            //end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
            end = config.lengthmz - 1 + nearfast(inp->dataMZ, (float)2 * inp->dataMZ[0] - point, config.lengthmz);
        }
        else {
            end = nearfast(inp->dataMZ, point, config.lengthmz);
        }
        endtab[i] = end;
        if (end - start > maxlength) { maxlength = end - start; }
        //printf("%d %d\n", start, end);
    }
    //printf("Max Length: %d\t", maxlength);
    return maxlength;
}


void cconv2fast(double* a, double* b, double* c, int length) {

    // We don't seem to have complex.h, unfortunately
    // fftw_complex is a double[2] of real (0) and imaginary (1)
    int complen = (length / 2) + 1;
    fftw_complex* A = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* B = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* C = fftw_malloc(complen * sizeof(fftw_complex));

    // Possible issue: create plan after we've created inputs (params)
    fftw_plan p1 = fftw_plan_dft_r2c_1d(length, a, A, FFTW_ESTIMATE);
    fftw_plan p2 = fftw_plan_dft_r2c_1d(length, b, B, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_execute(p2);

    // A * B = (ac - bd) + i(ad + bc)
    for (int j = 0; j < complen; j++) {
        C[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
        C[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
    }

    fftw_plan p3 = fftw_plan_dft_c2r_1d(length, C, c, FFTW_ESTIMATE);
    fftw_execute(p3); // Will destroy input array C, but that's okay

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_free(A);
    fftw_free(B);
    fftw_free(C);

}


void dd_deconv2(double* kernel_y, double* data_y, int length, double* output) {
    // Create flipped point spread function kernel_star
    double* kernel_star = malloc(length * sizeof(double));
    for (int i = 0; i < length; i++) {
        kernel_star[i] = kernel_y[length - i - 1];
    }
    // Create estimate for solution
    double* estimate = malloc(length * sizeof(double));
    for (int i = 0; i < length; i++) {
        estimate[i] = data_y[i];
    }
    // Allocate arrays for convolutions
    double* conv1 = malloc(length * sizeof(double));
    double* conv2 = malloc(length * sizeof(double));

    int complen = (length / 2) + 1;
    fftw_complex* kernel_ft = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* kernel_star_ft = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* estimate_ft = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* conv1_ft = fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* product_ft = fftw_malloc(complen * sizeof(fftw_complex));

    // FFTW_MEASURE takes a few seconds and overwrites input arrays, so we use ESTIMATE
    fftw_plan pk = fftw_plan_dft_r2c_1d(length, kernel_y, kernel_ft, FFTW_ESTIMATE); // for kernel_y
    fftw_plan pks = fftw_plan_dft_r2c_1d(length, kernel_star, kernel_star_ft, FFTW_ESTIMATE); // for kernel_star
    fftw_plan p1 = fftw_plan_dft_r2c_1d(length, estimate, estimate_ft, FFTW_ESTIMATE); // for estimate
    fftw_plan p2 = fftw_plan_dft_r2c_1d(length, conv1, conv1_ft, FFTW_ESTIMATE); // for conv1
    fftw_plan p1r = fftw_plan_dft_c2r_1d(length, product_ft, conv1, FFTW_ESTIMATE); // to conv1
    fftw_plan p2r = fftw_plan_dft_c2r_1d(length, product_ft, conv2, FFTW_ESTIMATE); // to conv2

    // We only need to do the transforms for kernel and kernel* once
    fftw_execute(pk);
    fftw_execute(pks);

    // Perform iterations
    int j = 0;
    double diff = 1.0;
    while (j < 50 && diff > 0.0001) { // Thresholds same as in Python
        // Find new estimate
        // cconv2fast(kernel_y, estimate, conv1, length);
        fftw_execute(p1);
        complex_mult(kernel_ft, estimate_ft, product_ft, complen);
        fftw_execute(p1r);
        for (int k = 0; k < length; k++) {
            if (conv1[k] != 0) {
                conv1[k] = data_y[k] / conv1[k];
            }
        }
        // cconv2fast(conv1, kernel_star, conv2, length);
        fftw_execute(p2);
        complex_mult(conv1_ft, kernel_star_ft, product_ft, complen);
        fftw_execute(p2r);
        for (int k = 0; k < length; k++) {
            conv2[k] = conv2[k] * estimate[k]; // Store new estimate in conv2
        }
        // Find how much the estimate changed
        double sum_diff = 0.0;
        double sum_est = 0.0;
        for (int k = 0; k < length; k++) {
            sum_diff += pow((estimate[k] - conv2[k]), 2.0);
            sum_est += estimate[k];
        }
        diff = sum_diff / sum_est;
        // Set new estimate as estimate
        for (int k = 0; k < length; k++) {
            estimate[k] = conv2[k];
        }
        j++;
    }

    // estimate now contains our deconvolution
    // Normalize and "return"
    double estimate_max = 0.0;
    for (int i = 0; i < length; i++) {
        if (estimate_max < estimate[i]) {
            estimate_max = estimate[i];
        }
    }
    for (int i = 0; i < length; i++) {
        output[i] = estimate[i] / estimate_max;
    }

    // Free arrays
    free(kernel_star);
    free(estimate);
    free(conv1);
    free(conv2);
    fftw_destroy_plan(pk);
    fftw_destroy_plan(pks);
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p1r);
    fftw_destroy_plan(p2r);
    fftw_free(kernel_ft);
    fftw_free(kernel_star_ft);
    fftw_free(estimate_ft);
    fftw_free(conv1_ft);
    fftw_free(product_ft);
}


// Convolution
