//
//  ErrorFunctions.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef ErrorFunctions_h
#define ErrorFunctions_h

#include <stdio.h>
#include "Config.h"
#include "Decon.h"
#include "FitFunctions.h"

__declspec(dllexport) float errfunspeedy(Config config, Decon decon, const char* barr, const float* dataInt, const int maxlength,
                   const int* isotopepos, const float* isotopeval, const int* starttab, const int* endtab, const float* mzdist, float* rsquared);


#endif /* ErrorFunctions_h */
