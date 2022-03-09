//
//  Transforms.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Transforms.h"

// Transforms //

void InterpolateTransform(const int maaxle, const int numz, const int lengthmz, const int* nztab, float* massaxis,
    const float adductmass, const float* dataMZ, float* massgrid, float* massaxisval, const float* blur)
{
    float startmzval = dataMZ[0];
    float endmzval = dataMZ[lengthmz - 1];
    //#pragma omp parallel for schedule(auto)
    for (int i = 0;i < maaxle;i++)
    {
        float val = 0;

        for (int j = 0;j < numz;j++)
        {
            float newval = 0;
            float mztest = (massaxis[i] + nztab[j] * adductmass) / (float)nztab[j];

            if (mztest > startmzval && mztest < endmzval)
            {
                int index = nearfast(dataMZ, mztest, lengthmz);
                int index2 = index;
                if (dataMZ[index] == mztest)
                {
                    newval = blur[index2D(numz, index, j)];
                    val += newval;
                    massgrid[index2D(numz, i, j)] = newval;
                }
                else
                {
                    if (dataMZ[index] > mztest && index > 1 && index < lengthmz - 1)
                    {
                        index2 = index;
                        index = index - 1;
                    }
                    else if (dataMZ[index] < mztest && index < lengthmz - 2 && index>0)
                    {
                        index2 = index + 1;
                    }
                    if (index2 > index && (dataMZ[index2] - dataMZ[index]) != 0)
                    {
                        float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                        float y0 = blur[index2D(numz, index - 1, j)];
                        float y1 = blur[index2D(numz, index, j)];
                        float y2 = blur[index2D(numz, index2, j)];
                        float y3 = blur[index2D(numz, index2 + 1, j)];
                        newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                        //newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
                        //newval=LinearInterpolate(y1,y2,mu);
                        val += newval;
                        massgrid[index2D(numz, i, j)] = newval;
                    }
                }
            }
        }
        massaxisval[i] = val;
    }
}

void SmartTransform(const int maaxle, const int numz, const int lengthmz, const int* nztab, float* massaxis,
    const float adductmass, const float* dataMZ, float* massgrid, float* massaxisval, const float* blur)
{
    float startmzval = dataMZ[0];
    float endmzval = dataMZ[lengthmz - 1];
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < maaxle; i++)
    {
        float val = 0;
        for (int j = 0; j < numz; j++)
        {
            float z = (float)nztab[j];
            float mtest = massaxis[i];
            float mztest = (mtest + z * adductmass) / z;

            float mzlower;
            float mlower;
            if (i > 0) {
                mlower = massaxis[i - 1];
                mzlower = (mlower + z * adductmass) / z;
                //mzlower = (mzlower + mztest) / 2;
            }
            else { mzlower = mztest; mlower = mtest; }

            float mzupper;
            float mupper;
            if (i < maaxle - 1) {
                mupper = massaxis[i + 1];
                mzupper = (mupper + z * adductmass) / z;
                //mzupper = (mzupper + mztest) / 2;
            }
            else { mzupper = mztest; mupper = mtest; }

            float newval = 0;
            if (mzupper > startmzval && mzlower < endmzval)
            {
                int index = nearfast(dataMZ, mztest, lengthmz);
                int index1 = nearfast(dataMZ, mzlower, lengthmz);
                int index2 = nearfast(dataMZ, mzupper, lengthmz);
                float imz = dataMZ[index];

                if (index2 - index1 < 5) {

                    if (imz == mztest)
                    {
                        newval = clip(blur[index2D(numz, index, j)], 0);
                        val += newval;
                        massgrid[index2D(numz, i, j)] = newval;

                    }
                    else
                    {
                        int edge = 0;
                        index2 = index;
                        if (imz > mztest)
                        {
                            index = index - 1;
                        }
                        else if (imz < mztest)
                        {
                            index2 = index + 1;
                        }

                        if (index < 1 || index2 >= lengthmz - 1)
                        {
                            edge = 1;
                        }
                        if (index < 0 || index2 >= lengthmz)
                        {
                            edge = 2;
                        }


                        if (index2 > index && (dataMZ[index2] - dataMZ[index]) != 0 && edge == 0)
                        {
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                            float y0 = blur[index2D(numz, index - 1, j)];
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = blur[index2D(numz, index2, j)];
                            float y3 = blur[index2D(numz, index2 + 1, j)];
                            newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                            //newval=clip(CRSplineInterpolate(y0,y1,y2,y3,mu), 0);
                            //newval=clip(LinearInterpolate(y1,y2,mu),0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("0\n");
                        }
                        else if (edge == 1 && (dataMZ[index2] - dataMZ[index]) != 0)
                        {
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = blur[index2D(numz, index2, j)];
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("1 %d %d %f %f %f\n", index, index2, mztest, imz, newval, massaxis[i]);
                        }
                        else if (edge == 2 && (dataMZ[index2] - dataMZ[index]) != 0)
                        {
                            float factor = 1;
                            if (index2 == 0) { index = 0; index2 = 1; }
                            if (index == lengthmz - 1) { index = lengthmz - 1; index2 = lengthmz - 2; factor = -1; }
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index] - dataMZ[index2]);
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = 0;
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("2\n");
                        }
                    }
                }

                else {
                    //printf("Integrating\n");
                    float num = 0;
                    for (int k = index1; k < index2 + 1; k++)
                    {
                        float kmz = dataMZ[k];
                        float km = (kmz - adductmass) * z;
                        float scale = 1;

                        if (mztest < kmz && km < mupper)
                        {
                            //scale = LinearInterpolatePosition(mzupper, mztest, kmz);
                            scale = LinearInterpolatePosition(mupper, mtest, km);
                        }
                        else if (kmz < mztest && km > mlower)
                        {
                            //scale = LinearInterpolatePosition(mzlower, mztest, kmz);
                            scale = LinearInterpolatePosition(mlower, mtest, km);
                        }
                        else if (kmz == mztest) { scale = 1; }
                        else { scale = 0; }

                        newval += scale * blur[index2D(numz, k, j)];
                        num += scale;

                    }
                    if (num != 0) { newval /= num; }
                    newval = clip(newval, 0);
                    val += newval;
                    massgrid[index2D(numz, i, j)] = newval;
                }
            }
        }
        massaxisval[i] = val;
    }
}



void IntegrateTransform(const int lengthmz, const int numz, const float* mtab, float massmax, float massmin,
    const int maaxle, float* massaxis, float* massaxisval, const float* blur, float* massgrid)
{
    for (int i = 0;i < lengthmz;i++)
    {
        for (int j = 0;j < numz;j++)
        {
            float testmass = mtab[index2D(numz, i, j)];
            if (testmass<massmax && testmass>massmin) {
                int index = nearfast(massaxis, testmass, maaxle);
                float newval = blur[index2D(numz, i, j)];
                if (massaxis[index] == testmass) {
                    massaxisval[index] += newval;
                    massgrid[index2D(numz, index, j)] += newval;
                }

                if (massaxis[index] < testmass && index < maaxle - 2)
                {
                    int index2 = index + 1;
                    float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
                    massaxisval[index] += (1.0 - interpos) * newval;
                    massgrid[index2D(numz, index, j)] += (1.0 - interpos) * newval;
                    massaxisval[index2] += (interpos)*newval;
                    massgrid[index2D(numz, index2, j)] += (interpos)*newval;
                }

                if (massaxis[index] > testmass && index > 0)
                {
                    int index2 = index - 1;
                    float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
                    massaxisval[index] += (1 - interpos) * newval;
                    massgrid[index2D(numz, index, j)] += (1 - interpos) * newval;
                    massaxisval[index2] += (interpos)*newval;
                    massgrid[index2D(numz, index2, j)] += (interpos)*newval;
                }
            }
        }
    }
}


 
// Transforms //
