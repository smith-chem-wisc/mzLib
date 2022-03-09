/*
* UD_analysis.h
*
*  Created on : 3 June 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//

#include "UD_dataproc.h"
#include "Config.h"

#ifndef ANALYSIS_HEADER
#define ANALYSIS_HEADER

void interpolate_merge(const float *massaxis, float *outint, const float *tempaxis, const float *tempint, const int mlen, const int templen);
void make_grid(int argc, char *argv[], Config config, const char *dtype, const char *out1, const char *out2, const char *out3);
int is_peak(const float *dataMZ, const float *dataInt, const int lengthmz, const float window, const float thresh, const int index);
int peak_detect(const float *dataMZ, const float *dataInt, const int lengthmz, const float window, const float thresh, float * peakx, float *peaky);
void peak_norm(float *peaky, int plen, int peaknorm);
float extract_height(Config config, const float peak, const float *xvals, const float *yvals, const int length);
float extract_localmax(Config config, const float peak, const float *xvals, const float *yvals, const int length);
float extract_localmax_position(Config config, const float peak, const float *xvals, const float *yvals, const int length);
float extract_integral(Config config, const float peak, const float *xvals, const float *yvals, const int length, const float thresh);
float extract_center_of_mass(Config config, const float peak, const float* xvals, const float* yvals, const int length, const float thresh); 
float extract_estimated_area(Config config, const float peak, const float* xvals, const float* yvals, const int length);
float extract_switch(Config config, const float peak, const float* xvals, const float* yvals, const int length);
void peak_extracts(Config config, const float *peakx, const char *dtype, int plen, int ultra);
void get_all_peaks(int argc, char *argv[], Config config);
void get_peaks(int argc, char *argv[], Config config, int ultra);


#endif
