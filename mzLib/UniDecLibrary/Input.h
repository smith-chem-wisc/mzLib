//
//  Input.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Input_h
#define Input_h

__declspec(dllexport) typedef struct Input Input;

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

__declspec(dllexport) typedef struct IsotopeStruct IsotopeStruct; 

struct IsotopeStruct {
    float minmid; 
    float minsig; 
    float maxmid; 
    float maxsig; 
    int isolength; 
};

__declspec(dllexport) Input SetupInputs(void);
__declspec(dllexport) void FreeInputs(Input inp);
__declspec(dllexport) Input SetupInputsSafe(Input* inp); 


#endif /* Input_h */
