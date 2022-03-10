//
//  Interpolation.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Interpolation.h"
#include <stdlib.h>
#include "Sorting.h"


//Perform a linear interpolation. I've left the code for other interpolation functions below but they don't seem to matter.
float LinearInterpolate(float y1, float y2, float mu)
{
    return(y1 * (1 - mu) + y2 * mu);
}

float LinearInterpolatePosition(float x1, float x2, float x)
{
    if (x2 - x1 == 0) { return 0; }
    return (x - x1) / (x2 - x1);
}

float CubicInterpolate(float y0, float y1, float y2, float y3, float mu)
{
    float a0, a1, a2, a3, mu2;
    mu2 = mu * mu;
    a0 = y3 - y2 - y0 + y1;
    a1 = y0 - y1 - a0;
    a2 = y2 - y0;
    a3 = y1;
    return(a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}

float CRSplineInterpolate(float y0, float y1, float y2, float y3, float mu)
{
    float a0, a1, a2, a3, mu2;
    mu2 = mu * mu;
    a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
    a1 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3;
    a2 = -0.5 * y0 + 0.5 * y2;
    a3 = y1;
    return(a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}

// Interpolation

double LinearInterpolateD(double y1, double y2, double mu)
{
    return(y1 * (1 - mu) + y2 * mu);
}

double LinearInterpolatePositionD(double x1, double x2, double x)
{
    if (x2 - x1 == 0) { return 0; }
    return (x - x1) / (x2 - x1);
}

// (Linear) Interpolates data in kernel to match sampling in data. Returns the new kernel length
int interpolate_dd(double* kernel_x, double* kernel_y, int kernellen, double* data_x,
    double* data_y, int datalen, double** kernel_x_new, double** kernel_y_new) {
    if (kernellen <= 1 || datalen <= 1) {
        return kernellen;
    }
    double diff = data_x[1] - data_x[0]; // kernel sampling needs to match this
    int newlen = ((kernel_x[kernellen - 1] - kernel_x[0]) / diff) + 1; // new number of points in kernel
    // these are the newly allocated arrays for the kernel
    int truelen;
    if (newlen > datalen) {
        truelen = newlen;
    } else {
        truelen = datalen;
    }
    *kernel_x_new = calloc(truelen, sizeof(double));
    *kernel_y_new = calloc(truelen, sizeof(double));

    double current_x = kernel_x[0];
    for (int i = 0; i < newlen; i++) {
        int nearest_index = nearfast_d(kernel_x, current_x, kernellen);
        if (kernel_x[nearest_index] == current_x) {
            (*kernel_y_new)[i] = kernel_y[nearest_index];
        } else if (kernel_x[nearest_index] < current_x) {
            if ((nearest_index + 1) < kernellen) {
                double mu = LinearInterpolatePositionD(kernel_x[nearest_index], kernel_x[nearest_index + 1], current_x);
                double y_val = LinearInterpolateD(kernel_y[nearest_index], kernel_y[nearest_index + 1], mu);
                (*kernel_y_new)[i] = y_val;
            } else { // this should never be the case
                (*kernel_y_new)[i] = kernel_y[nearest_index];
            }
        } else if (kernel_x[nearest_index] > current_x) {
            if (nearest_index > 0) {
                double mu = LinearInterpolatePositionD(kernel_x[nearest_index - 1], kernel_x[nearest_index], current_x);
                double y_val = LinearInterpolateD(kernel_y[nearest_index - 1], kernel_y[nearest_index], mu);
                (*kernel_y_new)[i] = y_val;
            } else { // this should also never happen
                (*kernel_y_new)[i] = kernel_y[nearest_index];
            }
        }
        (*kernel_x_new)[i] = current_x;
        current_x += diff;
    }

    return newlen;
}


// Interpolation

