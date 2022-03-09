//
//  Confidence.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Confidence.h"
// Confidence //


//Confidence score
float Confidence(float fit, float data)
{
    if (data == 0) { return 0; }
    return clip(1 - fabs(fit - data) / data, 0);
}


// Confidence //

