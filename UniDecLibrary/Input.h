//
//  Input.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Input_h
#define Input_h

typedef struct Input Input;

struct Input {
    float* dataMZ;
    float* dataInt;
    float* testmasses;
    int* nztab;
    float* mtab;
    char* barr;
    int* isotopepos;
    float* isotopeval;
    float isoparams[10];
};

Input SetupInputs(void);
void FreeInputs(Input inp);

#endif /* Input_h */
