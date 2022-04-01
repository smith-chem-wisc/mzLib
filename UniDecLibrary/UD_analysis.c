//
//  UD_analysis.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include <stdio.h>
#include <stdlib.h>
#include "MathUtilities.h"
#include "Config.h"
#include "Decon.h"
#include "Input.h"
#include "Sorting.h"
#include "Interpolation.h"
#include "UD_score.h"

void interpolate_merge(const float *massaxis, float *outint, const float *tempaxis, const float *tempint, const int mlen, const int templen)
{
    float start = tempaxis[0];
    float end = tempaxis[templen - 1];

    for (int i = 0; i < mlen; i++)
    {
        float pos = massaxis[i];
        if (pos >= start && pos <= end)
        {
            int index = nearfast(tempaxis, pos, templen);
            int index2 = index;
            if (tempaxis[index] == pos)
            {
                outint[i] = tempint[index];
            }
            else
            {
                if (tempaxis[index]>pos&&index>1 && index<templen - 1)
                {
                    index2 = index;
                    index = index - 1;
                }
                else if (tempaxis[index]<pos&&index<templen - 2 && index>0)
                {
                    index2 = index + 1;
                }
                if (index2>index && (tempaxis[index2] - tempaxis[index]) != 0)
                {
                    float mu = (pos - tempaxis[index]) / (tempaxis[index2] - tempaxis[index]);
                    float y0 = tempint[index - 1];
                    float y1 = tempint[index];
                    float y2 = tempint[index2];
                    float y3 = tempint[index2 + 1];
                    float newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                    //newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
                    //newval=LinearInterpolate(y1,y2,mu);
                    outint[i] = newval;
                }
            }
        }
    }

}


int is_peak(const float* dataMZ, const float* dataInt, const int lengthmz, const float window, const float thresh, const int index)
{
    float xval = dataMZ[index];
    float yval = dataInt[index];
    //Check if below threshold
    if (yval < thresh) { return 0; }
    //Check if local max
    for (int i = 0; i < lengthmz; i++)
    {
        float temp = dataMZ[i];
        if (fabs(temp - xval) <= window)
        {
            float tempy = dataInt[i];
            //printf("%f %f %f %f\n", temp, tempy, xval, yval);
            if (tempy > yval)
            {
                return 0;
            }
            if (tempy == yval && i < index)
            {
                return 0;
            }
        }
    }

    return 1;
}





int peak_detect(const float *dataMZ, const float *dataInt, const int lengthmz, const float window, const float thresh, float * peakx, float *peaky)
{
    //printf("Detecting Peaks %d %f %f\n", lengthmz, window, thresh);
    int plen = 0;
    float max = Max(dataInt, lengthmz);
    for (int i = 0; i < lengthmz; i++)
    {
        if (is_peak(dataMZ, dataInt, lengthmz, window, thresh*max, i) == 1)
        {
            //printf("Peak %d: %f %f\n", plen, dataMZ[i], dataInt[i]);//
            peakx[plen] = dataMZ[i];
            peaky[plen] = dataInt[i];
            plen++;
        }
    }
    return plen;
}


__declspec(dllexport) void peak_norm(float *peaky, int plen, int peaknorm)
{
    float norm = 0;
    if (peaknorm==1){norm= Max(peaky, plen);}
    if (peaknorm==2){norm= Sum(peaky, plen);}
    if (norm != 0)
    {
        for (int i = 0; i < plen; i++)
        {
            peaky[i] /= norm;
        }
    }
}

float extract_height(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
    if (peak < xvals[0]) { return 0; }
    if (peak > xvals[length - 1]) { return 0; }
    int pos = nearfast(xvals, peak, length);
    return yvals[pos];
}

float extract_localmax(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
    if (peak < xvals[0]) { return 0; }
    if (peak > xvals[length - 1]) { return 0; }

    int pos1 = nearfast(xvals, peak - config.exwindow, length);
    int pos2 = nearfast(xvals, peak + config.exwindow, length);

    float localmax = 0;
    for (int i = pos1; i <= pos2; i++)
    {
        if (yvals[i] > localmax) { localmax = yvals[i]; }
    }

    return localmax;
}

float extract_localmax_position(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
    if (peak < xvals[0]) { return 0; }
    if (peak > xvals[length - 1]) { return 0; }

    int pos1 = nearfast(xvals, peak - config.exwindow, length);
    int pos2 = nearfast(xvals, peak + config.exwindow, length);

    float localmax = 0;
    int localmaxpos = 0;
    for (int i = pos1; i <= pos2; i++)
    {
        if (yvals[i] > localmax) { localmax = yvals[i]; localmaxpos = i; }
    }

    return xvals[localmaxpos];
}

float extract_integral(Config config, const float peak, const float *xvals, const float *yvals, const int length, const float thresh)
{
    if (peak < xvals[0]) { return 0; }
    if (peak > xvals[length - 1]) { return 0; }

    float thresh2 = 0;
    if (thresh > 0) {
        float max = Max(yvals, length);
        thresh2 = thresh * max;
    }
    //printf("thresh %f\n", thresh2);

    int pos1 = nearfast(xvals, peak - config.exwindow, length);
    int pos2 = nearfast(xvals, peak + config.exwindow, length);

    float integral = 0;
    for (int i = pos1+1; i <= pos2; i++)
    {
        float b = xvals[i];
        float a = xvals[i - 1];
        float fb = yvals[i];
        float fa = yvals[i - 1];
        if (fa > thresh2 && fb > thresh2) {
            integral += (b - a) * ((fa + fb) / 2.0);
        }
    }

    return integral;
}


float extract_center_of_mass(Config config, const float peak, const float *xvals, const float *yvals, const int length, const float thresh)
{
    if (peak < xvals[0]) { return 0; }
    if (peak > xvals[length - 1]) { return 0; }
        
    float thresh2 = 0;
    if (thresh > 0) {
        float max = Max(yvals, length);
        thresh2 = thresh * max;
    }
    //printf("thresh %f\n", thresh2);

    int pos1 = nearfast(xvals, peak - config.exwindow, length);
    int pos2 = nearfast(xvals, peak + config.exwindow, length);
    float sum = 0;
    float sum_masses = 0;
    for (int i = pos1; i <= pos2; i++)
    {
        float x = xvals[i];
        float y = yvals[i];
        if (y > thresh2)
        {
            sum_masses += y;
            sum += x*y;
        }
    }
    if (sum_masses > 0) { sum /= sum_masses; }

    return sum;
}


float extract_estimated_area(Config config, const float peak, const float* xvals, const float* yvals, const int length)
{
    int pos1 = nearfast(xvals, peak - config.exwindow, length);
    // printf("Position of left bound is %d\n", pos1);
    int pos2 = nearfast(xvals, peak + config.exwindow, length);
    // printf("Position of right bound is %d\n", pos2);
    int mlen = pos2 - pos1 + 1;
    const float* xwin = xvals + pos1;
    // printf("xvals is %p, xwin is %p\n", xvals, xwin);
    const float* ywin = yvals + pos1;

    int index = nearfast(xwin, peak, mlen);
    float height = ywin[index];
    float fwhm = single_fwhm(config, mlen, xwin, ywin, peak, index, height);
    
    float pi = 3.14159265358979323846;
    float gauss_coeff = sqrt(pi / log(2.0)) / 2.0;
    float adjusted_coeff = ((0.5 * gauss_coeff) + (pi / 4.0));
    float area = 0;
    if (config.psfun == 0) { // Gaussian
        area = height * fwhm * gauss_coeff;
    }
    else if (config.psfun == 1) { // Lorentzian
        area = height * fwhm * pi / 2.0;
    }
    else if (config.psfun == 2) { // Split G/L
        area = height * fwhm * adjusted_coeff;
    }

    return area;
}


float extract_switch(Config config, const float peak, const float* xvals, const float* yvals, const int length)
{
    float output = 0;
    if (config.exwindow == 0) { config.exchoice = 0; }
    
    float thresh = config.exthresh / 100;

    switch (config.exchoice)
    {
    case 0:
        //printf("Extracting Height\n");
        output = extract_height(config, peak, xvals, yvals, length);
        break;
    case 1:
        //printf("Extracting Local Max. Window: %f\n", config.exwindow);
        output = extract_localmax(config, peak, xvals, yvals, length);
        break;
    case 2:
        //printf("Extracting Integral. Window: %f\n", config.exwindow);
        output = extract_integral(config, peak, xvals, yvals, length, thresh);
        break;
    case 3:
        //printf("Extracting Center of Mass 0. Window: %f\n", config.exwindow);
        output = extract_center_of_mass(config, peak, xvals, yvals, length, thresh);
        break;
    case 4:
        //printf("Extracting Local Max. Window: %f\n", config.exwindow);
        output = extract_localmax_position(config, peak, xvals, yvals, length);
        break;
    /*
    case 5:
        //printf("Extracting Center of Mass 50. Window: %f\n", config.exwindow);
        output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.5*max);
        break;
    case 6:
        //printf("Extracting Center of Mass 10. Window: %f\n", config.exwindow);
        output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.1*max);
        break;
    case 7:
        //printf("Extracting Integral. Window: %f\n", config.exwindow);
        output = extract_integral(config, peak, xvals, yvals, length, 0.5 * max);
        break;
    case 8:
        //printf("Extracting Integral. Window: %f\n", config.exwindow);
        output = extract_integral(config, peak, xvals, yvals, length, 0.1 * max);
        break;
    case 9:
        //printf("Extracting Integral. Window: %f\n", config.exwindow);
        output = extract_integral(config, peak, xvals, yvals, length, 0.05 * max);
           break;
    case 10:
        //printf("Extracting Integral. Window: %f\n", config.exwindow);
        output = extract_integral(config, peak, xvals, yvals, length, 0.025 * max);
        break;*/
    case 5: // 11
        // printf("Extracting Estimated Area. Window: %f\n", config.exwindow);
        output = extract_estimated_area(config, peak, xvals, yvals, length);
        break;
    default:
        printf("Invalid Extraction Choice: %d\n", config.exchoice);
        output = 0;
    }
    return output;
}

