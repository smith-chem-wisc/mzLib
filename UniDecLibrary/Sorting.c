//
//  Sorting.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Sorting.h"

// Sorting functions

//Slow nearest unsorted.
int nearunsorted(float* testmasses, float point, int lengthtest)
{
    float minval = fabs(point - testmasses[0]);
    float val = testmasses[0];
    float difftest;
    int pos = 0;
    for (int i = 1; i < lengthtest; i++)
    {
        difftest = fabs(point - testmasses[i]);
        if (difftest < minval)
        {
            minval = difftest;
            val = testmasses[i];
            pos = i;
        }
    }
    return pos;
}

//Slow way to test if two points are within a cutoff distance of each other. Works on unsorted list.
int neartest(float* testmasses, float point, int lengthtest, float cutoff)
{
    float minval = fabs(point - testmasses[0]);
    float val = testmasses[0];
    for (int i = 0; i < lengthtest; i++)
    {
        float difftest = fabs(point - testmasses[i]);
        if (difftest < minval)
        {
            minval = difftest;
            val = testmasses[i];
        }
    }
    int test = 0;
    if (fabs(point - val) < cutoff)
    {
        test = 1;
    }
    return test;
}

//Fast way of finding the nearest data point in an ordered list.
int nearfast(const float* dataMZ, const float point, const int numdat)
{
    int start = 0;
    int length = numdat - 1;
    int end = 0;
    int diff = length - start;
    while (diff > 1)
    {
        if (point < dataMZ[start + (length - start) / 2])
        {
            length = start + (length - start) / 2;
        }
        else if (point == dataMZ[start + (length - start) / 2])
        {
            end = start + (length - start) / 2;
            length = start + (length - start) / 2;
            start = start + (length - start) / 2;
            return end;
        }
        else if (point > dataMZ[start + (length - start) / 2])
        {
            start = start + (length - start) / 2;
        }
        diff = length - start;
    }
    if (fabs(point - dataMZ[start]) >= fabs(point - dataMZ[length]))
    {
        end = length;
    }
    else
    {
        end = start;
    }
    return end;
}

int nearfast_test(const float* dataMZ, const float point, const int numdat, float cutoff)
{
    int index = nearfast(dataMZ, point, numdat);
    float val = dataMZ[index];
    if (fabs(val - point) < cutoff)
    {
        return index;
    }
    else { return -1; }
}


// Fast way of finding the nearest data point in an ordered list of doubles.
int nearfast_d(const double* dataMZ, const double point, const int numdat)
{
    int start = 0;
    int length = numdat - 1;
    int end = 0;
    int diff = length - start;
    while (diff > 1) {
        if (point < dataMZ[start + (length - start) / 2]) {
            length = start + (length - start) / 2;
        } else if (point == dataMZ[start + (length - start) / 2]) {
            end = start + (length - start) / 2;
            length = start + (length - start) / 2;
            start = start + (length - start) / 2;
            return end;
        } else if (point > dataMZ[start + (length - start) / 2]) {
            start = start + (length - start) / 2;
        }
        diff = length - start;
    }
    if (fabs(point - dataMZ[start]) >= fabs(point - dataMZ[length])) {
        end = length;
    } else {
        end = start;
    }
    return end;
}


// Sorting
