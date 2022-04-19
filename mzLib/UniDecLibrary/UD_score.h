/*
* UD_score.h
*
*  Created on : 20 August 2019
* Author : Michael.Marty
*/

//
// 
// Copyright 2019 University of Arizona
//
//

#include "UD_analysis.h"
#include "Config.h"
#include "Input.h"
#include "Decon.h"
//#include "UniDec.h"
//#include "UD_peak_width.h"
void get_fwhms(Config config, const int plen, const int mlen, const float* massaxis, const float* masssum, const float* peakx, float* fwhmlow, float* fwhmhigh, float* badfwhm);
float single_fwhm(Config config, const int mlen, const float* massaxis, const float* masssum, const float peak, int index, float max);
float uscore(Config config, const float* dataMZ, const float* dataInt, const float* mzgrid, const int* nztab,
             const float mlow, const float mhigh, const float peak); 
float mscore(Config config, const int mlen, const float* massaxis, const float* masssum, const float* massgrid, const float mlow, const float mhigh, const float peak);
float csscore(Config config, const int mlen, const float* massaxis, const float* masssum, const float* massgrid, const float mlow, const float mhigh, const float peak);
float find_minimum(Config config, const int mlen, const float* massaxis, const float* masssum, const float lowpt, const float highpt);
float score_minimum(float height, float min);
float fscore(Config config, const int plen, const int mlen, const float* massaxis, const float* masssum, const float *peakx, const float height,
             const float mlow, const float mhigh, const float peak, const int badfwhm);
float score_from_peaks(const int plen, const float *peakx, const float *peaky, float *dscores, const Config config, Decon *decon, const Input inp, const float threshold);
float performScoring(Config config, Decon *decon, Input inp, const float threshold);
//int ReadDecon(Config* config, const Input inp, Decon* decon);
//void get_scan_scores(int argc, char* argv[], Config config);

