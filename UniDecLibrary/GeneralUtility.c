//
//  GeneralUtility.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "GeneralUtility.h"

// General Utility //


void deepcopy(float* out, const float* in, int len)
{
    for (int i = 0; i < len; i++)
    {
        out[i] = in[i];
    }
}

