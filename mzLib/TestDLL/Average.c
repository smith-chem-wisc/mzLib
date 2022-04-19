#include "Average.h"
float Average(const int length, const float* xarray)
{
    float temp1 = 0;
    float temp2 = (float)length;
    for (int i = 0; i < length; i++)
    {
        temp1 += xarray[i];
    }
    if (temp2 == 0) { return 0; }
    return temp1 / temp2;
}