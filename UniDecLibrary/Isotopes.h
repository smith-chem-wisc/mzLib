//
//  Isotopes.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Isotopes_h
#define Isotopes_h

#include <stdio.h>
#include "Config.h"
#include "Input.h"
#include "ArrayIndexing.h"

void setup_and_make_isotopes(Config* config, Input* inp);
void test_isotopes(float mass, float* isoparams);
void isotope_dist(float mass, int isolength, int* isoindex, float* isovals, float* isoparams);
void make_isotopes(float* isoparams, int* isotopepos, float* isotopeval, float* mtab, int* ztab, char* barr, float* dataMZ, int lengthmz, int numz);
int setup_isotopes(float* isoparams, int* isotopepos, float* isotopeval, float* mtab, int* ztab, char* barr, float* dataMZ, int lengthmz, int numz);
float isotopebeta(float mass, float* isoparams);
float isotopealpha(float mass, float* isoparams);
float isotopesig(float mass, float* isoparams);
float isotopemid(float mass, float* isoparams); 
void monotopic_to_average(const int lengthmz, const int numz, float* blur, const char* barr, int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval);
#endif /* Isotopes_h */
