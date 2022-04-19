//
//  Interpolation.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Interpolation_h
#define Interpolation_h

#include <stdio.h>

float LinearInterpolate(float y1, float y2, float mu);
float LinearInterpolatePosition(float x1, float x2, float x);
float CubicInterpolate(float y0, float y1, float y2, float y3, float mu);
float CRSplineInterpolate(float y0, float y1, float y2, float y3, float mu);
double LinearInterpolateD(double y1, double y2, double mu);
double LinearInterpolatePositionD(double x1, double x2, double x);
int interpolate_dd(double* kernel_x, double* kernel_y, int kernellen, double* data_x,
                   double* data_y, int datalen, double** kernel_x_new, double** kernel_y_new);

#endif /* Interpolation_h */
