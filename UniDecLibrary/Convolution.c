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
        /*
        for (i = 0; i<lengthmz; i++)
        {
            deltas[lengthmz - 1 - i] = denom[i];
        }
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
        for (i = 0; i<lengthmz; i++)
        {
            deltas[lengthmz - 1 - i] = denom[i];
        }
        memcpy(denom, deltas, sizeof(float)*lengthmz);*/
    }

    //printf("6\n");
    //Multiply Ratio by prior
    apply_ratios(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, denom, blur2);

    //printf("7\n");
    if (aggressiveflag == 1)
    {
        //memcpy(deltas, denom, sizeof(float)*lengthmz);
        blur_baseline(denom, lengthmz, dataMZ, fabs(mzsig), 0, filterwidth);
        //blur_baseline(denom, lengthmz, 10);
        //blur_noise(deltas, lengthmz);
#pragma omp parallel for private(i), schedule(auto)
        for (i = 0; i < lengthmz; i++)
        {
            baseline[i] = baseline[i] * (denom[i]);
            //noise[i] = noise[i]*(deltas[i]);
        }
    }
    //printf("8\n");
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
            start = 0 - nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
        }
        else {
            start = nearfast(inp->dataMZ, point, config.lengthmz);
        }
        starttab[i] = start;

        point = inp->dataMZ[i] + threshold;
        if (point > inp->dataMZ[config.lengthmz - 1] && config.speedyflag == 0) {
            //end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
            end = config.lengthmz - 1 + nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
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


// Gives convolution of functions a and b. Unused.
// Use cconv2fast instead
void cconv2(double* a, double* b, double* c, int length) {
    double** A = malloc(length * sizeof(double*));    // (a + bi)
    double** B = malloc(length * sizeof(double*));    // (c + di)
    double** C = malloc(length * sizeof(double*));
    for (int i = 0; i < length; i++) {
        A[i] = calloc(2, sizeof(double));
        B[i] = calloc(2, sizeof(double));
        C[i] = calloc(2, sizeof(double));
    }

    discretefouriertransform(a, A, length);
    discretefouriertransform(b, B, length);
    // A * B = (ac - bd) + i(ad + bc)
    for (int j = 0; j < length; j++) {
        C[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
        C[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
    }
    inversefouriertransform(C, c, length);
    for (int k = 0; k < length; k++) {
        if (c[k] < 0) {
            c[k] = c[k] * -1.0;
        }
        free(A[k]);
        free(B[k]);
        free(C[k]);
    }
    free(A);
    free(B);
    free(C);
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

void dd_deconv(double* kernel_y, double* data_y, int length, double* output) {

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

    // Perform iterations
    int j = 0;
    double diff = 1.0;
    while (j < 50 && diff > 0.0001) { // Thresholds same as in Python
        // Find new estimate
        cconv2fast(kernel_y, estimate, conv1, length);
        for (int k = 0; k < length; k++) {
            if (conv1[k] != 0) {
                conv1[k] = data_y[k] / conv1[k];
            }
        }
        cconv2fast(conv1, kernel_star, conv2, length);
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

}


void DoubleDecon(const Config* config, Decon* decon) {

    // Get larger length
    int true_length;
    int kernel_length = getfilelength(config->kernel);
    int data_length = decon->mlen;
    if (kernel_length > data_length) {
        true_length = kernel_length;
    } else {
        true_length = data_length;
    }

    // Read in kernel file
    double* kernel_x_init = calloc(true_length, sizeof(double));
    double* kernel_y_init = calloc(true_length, sizeof(double));
    readkernel(config->kernel, kernel_length, kernel_x_init, kernel_y_init);
    // printf("Kernel file read.\n");

    // Enforce the same sampling on the kernel
    if (kernel_length > 1 && data_length > 1) {
        double diff = decon->massaxis[1] - decon->massaxis[0]; // kernel sampling needs to match this
        double kdiff = kernel_x_init[1] - kernel_x_init[0]; // the original kernel sampling
        if (diff != kdiff) {
            int newlen = ((kernel_x_init[kernel_length - 1] - kernel_x_init[0]) / diff) + 1;
            if (newlen > data_length) {
                true_length = newlen;
            } else {
                true_length = data_length;
            }
        }
    }

    // Read in data (i.e. copy from decon struct)
    double* data_x = calloc(true_length, sizeof(double));
    double* data_y = calloc(true_length, sizeof(double));
    double max_data_y = 0.0;
    for (int i = 0; i < data_length; i++) {
        data_x[i] = decon->massaxis[i];
        data_y[i] = decon->massaxisval[i];
        if (data_y[i] > max_data_y) {
            max_data_y = data_y[i];
        }
    }
    for (int i = 0; i < data_length; i++) {    // Normalize
        data_y[i] = data_y[i] / max_data_y;
    }
    // printf("Data file copied and normalized.\n");

    // Integrate or interpolate if necessary...
    double* kernel_x = NULL;
    double* kernel_y = NULL;
    int kernel_length2 = true_length;
    // printf("true length is %d\n", kernel_length2);
    if (kernel_length > 1 && data_length > 1 &&
        (data_x[1] - data_x[0]) > (kernel_x_init[1] - kernel_x_init[0])) {
        // printf("Integrating\n");
        kernel_length2 = integrate_dd(kernel_x_init, kernel_y_init, kernel_length, data_x, data_y,
            true_length, &kernel_x, &kernel_y);
    } else if (kernel_length > 1 && data_length > 1 &&
        (data_x[1] - data_x[0]) < (kernel_x_init[1] - kernel_x_init[0])) {
        // printf("Interpolating\n");
        kernel_length2 = interpolate_dd(kernel_x_init, kernel_y_init, kernel_length, data_x, data_y,
            true_length, &kernel_x, &kernel_y);
    } else {
        // printf("Sampling is OK\n");
        kernel_x = kernel_x_init;
        kernel_y = kernel_y_init;
    }
    // ...and find max and normalize
    double max_kernel_y = 0.0;
    int max_kernel_i = 0;
    for (int i = 0; i < kernel_length2; i++) {
        // printf("kernel y for index %d is %f\n", i, kernel_y[i]);
        if (kernel_y[i] > max_kernel_y) {
            max_kernel_y = kernel_y[i];
            max_kernel_i = i;
        }
    }
    // printf("max is %f, out of length %d\n", max_kernel_y, kernel_length2);
    for (int i = 0; i < kernel_length2; i++) {    // Normalize
        kernel_y[i] = kernel_y[i] / max_kernel_y;
        // printf("%f\n", kernel_y[i]);
    }
    // printf("Kernel file normalized.\n");

    // Pad x-axis for the shorter one (is padding kernel even necessary???)
    if (data_length < true_length) { // Pad data_x
        double diff = data_x[1] - data_x[0];
        double last = data_x[data_length - 1];
        for (int i = data_length; i < true_length; i++) {
            data_x[i] = last + diff;
            last = data_x[i];
        }
    }
    else if (kernel_length < true_length) { // Pad kernel_x
        double diff = kernel_x[1] - kernel_x[0];
        double last = kernel_x[kernel_length - 1];
        for (int i = kernel_length; i < true_length; i++) {
            kernel_x[i] = last + diff;
            last = kernel_x[i];
        }
    }
    // printf("Data padded.\n");

    // Prepare kernel
    double* real_kernel_y = calloc(true_length, sizeof(double));
    int part1_length = true_length - max_kernel_i;
    for (int i = 0; i < part1_length; i++) {
        real_kernel_y[i] = kernel_y[max_kernel_i + i];
    }
    for (int i = 0; i < max_kernel_i; i++) {
        real_kernel_y[part1_length + i] = kernel_y[i];
    }
    // printf("Kernel file prepared.\n");

    // Run Richardson-Lucy deconvolution
    double* doubledec = calloc(true_length, sizeof(double));
    // printf("Running dd_deconv2\n");
    dd_deconv2(real_kernel_y, data_y, true_length, doubledec);

    // printf("Copying results to Decon struct.\n");
    int lb = nearfast_d(data_x, config->masslb, true_length);
    if (data_x[lb] < config->masslb) lb++;
    int ub = nearfast_d(data_x, config->massub, true_length);
    if (data_x[ub] > config->massub) lb--;
    int write_length = ub - lb + 1;
    if (write_length > decon->mlen) {
        // printf("Warning: new length exceeds previous mlen.\n");
        free(decon->massaxis);
        free(decon->massaxisval);
        decon->massaxis = calloc(write_length, sizeof(float));
        decon->massaxisval = calloc(write_length, sizeof(float));
    }
    // Copy results to the Decon struct
    for (int i = 0; i < write_length; i++) {
        decon->massaxis[i] = data_x[i + lb];
        decon->massaxisval[i] = doubledec[i + lb];
    }
    decon->mlen = write_length;
    // printf("Results copied to Decon.\n");

    // Free memory
    if (kernel_x != kernel_x_init) {
        free(kernel_x);
        free(kernel_y);
    }
    free(kernel_x_init);
    free(kernel_y_init);
    free(real_kernel_y);
    free(data_x);
    free(data_y);
    free(doubledec);

}


// Convolution
