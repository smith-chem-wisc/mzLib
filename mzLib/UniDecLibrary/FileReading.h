//
//  FileReading.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef FileReading_h
#define FileReading_h

#include <stdio.h>

void readfile(char* infile, int lengthmz, float* dataMZ, float* dataInt);
void readfile3(char* infile, int lengthmz, float* array1, float* array2, float* array3);
void readfile3bin(char* infile, int lengthmz, float* array1, float* array2, float* array3);
void readmfile(char* infile, int mfilelen, float* testmasses);
int getfilelength(const char* infile);
int getfilelengthbin(const char* infile, const int size, const int width);

#endif /* FileReading_h */
