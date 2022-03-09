//
//  Integration.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Integration.h"

// Integration //

// Integrates data in kernel to match sampling in data. Returns the new kernel length
int integrate_dd(double* kernel_x, double* kernel_y, int kernellen, double* data_x,
    double* data_y, int datalen, double** kernel_x_new, double** kernel_y_new) {
    if (kernellen <= 1 || datalen <= 1) {
        return kernellen;
    }
    double diff = data_x[1] - data_x[0]; // kernel sampling needs to match this
    double kdiff = kernel_x[1] - kernel_x[0]; // the original kernel sampling
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
    double current_xl = kernel_x[0];
    double current_xr = kernel_x[0] + (diff / 2);
    int current_index = 0;
    for (int i = 0; i < newlen; i++) {
        double y_val = 0;
        for (int j = current_index; j < kernellen; j++) {
            // For the first value, add area to the left of the point
            if (j == current_index && j != 0 && kernel_x[j] >= current_xl && kernel_x[j] < current_xr) {
                double left_mu = LinearInterpolatePositionD(kernel_x[j - 1], kernel_x[j], current_xl);
                double left_y = LinearInterpolateD(kernel_y[j - 1], kernel_y[j], left_mu);
                y_val += (left_y + kernel_y[j]) * (kernel_x[j] - current_xl) / 2.0;
            }
            // Next, add the area to the right of the point (it's either to the next point, to the
            // boundary, or we're at the last point)
            if (kernel_x[j] >= current_xl && kernel_x[j] < current_xr && (j + 1) < kernellen &&
                kernel_x[j + 1] < current_xr) {
                y_val += (kernel_y[j] + kernel_y[j + 1]) * kdiff / 2.0;
            } else if (kernel_x[j] >= current_xl && kernel_x[j] < current_xr &&
                (j + 1) < kernellen && kernel_x[j + 1] >= current_xr) {
                double right_mu = LinearInterpolatePositionD(kernel_x[j], kernel_x[j + 1], current_xr);
                double right_y = LinearInterpolateD(kernel_y[j], kernel_y[j + 1], right_mu);
                y_val += (kernel_y[j] + right_y) * (current_xr - kernel_x[j]) / 2.0;
            } else if (kernel_x[j] >= current_xr || (j + 1) >= kernellen) {
                current_index = j;
                break;
            }
        }
        (*kernel_x_new)[i] = current_x;
        (*kernel_y_new)[i] = y_val; // Should probably divide by diff. CHECK! --actually, I don't think it matters
        current_x += diff;
        current_xl = current_xr;
        current_xr += diff;
    }

    return newlen;
}
// Integration //
