#pragma once
/*
* UD_dataproc.h
*
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//
#include <stdlib.h>
#include "Config.h"

#ifndef DATAPROC_HEADER
#define DATAPROC_HEADER

int pool1d(float *oldmz, float *oldint, const int oldlen, const int mzbins);
int chop1d(float *oldmz, float *oldint, int oldlen, const float min, const float max); 
int remove_duplicates(float *oldmz, float *oldint, int oldlen);
int remove_middle_zeros(float *oldmz, float *oldint, int oldlen);
void norm1d(float *dataInt, const int lengthmz);
void background_subtract(float *dataInt, const float bsub, const int lengthmz);
inline int compare(const void* a, const void* b);
int data_reduction(float* oldmz, float* oldint, int oldlen, const float redper);
void process_data(int argc, char *argv[], Config config);

#endif
