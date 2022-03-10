//
//  Normalization.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Normalization.h"
#include "MathUtilities.h"
// Normalization //

void simp_norm(const int length, float* data)
{
    float norm = Max(data, length);
    //printf("\nMax: %f %d\n", norm, length);
    if (norm > 0) {
        for (int i = 0; i < length; i++)
        {
            data[i] = data[i] / norm;
        }
    }
    return;
}

void simp_norm_sum(const int length, float* data)
{
    float norm = Sum(data, length);
    //printf("\nMax: %f %d\n", norm, length);
    if (norm > 0) {
        for (int i = 0; i < length; i++)
        {
            data[i] = data[i] / norm;
        }
    }
    return;
}
