//
//  Unused.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Unused_h
#define Unused_h

#include <stdio.h>
void IntPrint(const int* array, const int length);
void floatPrint(const float* array, const int length);
void textvectorprint(float* arr, int length);
void discretefouriertransform(double* input, double** output, int length);
void inversefouriertransform(double** input, double* output, int length);
void WriteDecon(const Config config, const Decon* decon, const Input* inp);
void ReadInputs(int argc, char* argv[], Config* config, Input* inp);
void WritePeaks(const Config config, const Decon* decon);
void SetLimits(const Config config, Input* inp);
void monotopic_to_average(const int lengthmz, const int numz, float* blur, const char* barr, int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval);
Config PostImport(Config config);
void ManualAssign(float* dataMZ, char* barr, int* nztab, Config config);

#endif /* Unused_h */
