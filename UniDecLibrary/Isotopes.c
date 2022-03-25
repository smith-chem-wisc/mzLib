//
//  Isotopes.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Isotopes.h"
#include "Sorting.h"



// Isotopes //


float isotopemid(float mass, float* isoparams)
{
    float a, b, c;
    a = isoparams[4];
    b = isoparams[5];
    c = isoparams[6];
    return a + b * pow(mass, c);
}

float isotopesig(float mass, float* isoparams)
{
    float a, b, c;
    a = isoparams[7];
    b = isoparams[8];
    c = isoparams[9];
    return a + b * pow(mass, c);
}

float isotopealpha(float mass, float* isoparams)
{
    float a, b;
    a = isoparams[0];
    b = isoparams[1];
    return a * exp(-mass * b);
}

float isotopebeta(float mass, float* isoparams)
{
    float a, b;
    a = isoparams[2];
    b = isoparams[3];
    return a * exp(-mass * b);
}

int setup_isotopes(float* isoparams, int* isotopepos, float* isotopeval, float* mtab, int* ztab, char* barr, float* dataMZ, int lengthmz, int numz)
{
    float minmass = 100000000;
    float maxmass = 1;
    int i;
    for (i = 0;i < lengthmz * numz;i++)
    {
        if (barr[i] == 1)
        {
            float mass = mtab[i];
            if (mass < minmass) { minmass = mass; }
            if (mass > maxmass) { maxmass = mass; }
        }
    }

    float minmid = isotopemid(minmass, isoparams);
    float minsig = isotopesig(minmass, isoparams);
    float maxmid = isotopemid(maxmass, isoparams);
    float maxsig = isotopesig(maxmass, isoparams);

    int isostart = (int)(minmid - 4 * minsig);
    int isoend = (int)(maxmid + 4 * maxsig);
    if (isostart < 0) { isostart = 0; }
    if (isoend < 4) { isoend = 4; }
    int isolength = isoend - isostart;
    return isolength;
}

// use a pointer to pass minmid, maxmid, maxsig. 
void make_isotopes(float* isoparams, int* isotopepos, float* isotopeval, float* mtab, int* ztab, char* barr,
    float* dataMZ, int lengthmz, int numz, float minmid, float maxmid, float maxsig)
{
    float massdiff = 1.0026;
    int i, j, k;

    int isostart = 0;
    int isoend = (int)(maxmid + 4 * maxsig);
    //if (isostart<0){ isostart = 0; }
    if (isoend < 4) { isoend = 4; }
    
    int isolength = isoend - isostart;
    
    float* isorange = NULL;
    int* isoindex = NULL;
    isorange = calloc(isolength, sizeof(float));
    isoindex = calloc(isolength, sizeof(int));
    // looping through and assigning isorange and isoindex. 
    for (i = 0; i < isolength; i++)
    {
        // should probably fix this null dereferencing. 
        isorange[i] = (isostart + i) * massdiff;
        isoindex[i] = (isostart + i);
    }
    // iterate and get rid of the non-zero values, then iterate over that. saves 
    // an if for each iteration. 
#pragma omp parallel for private (i,j,k), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        for (j = 0; j < numz; j++)
        {
            if (barr[index2D(numz, i, j)] == 1)
            {
                float mz = dataMZ[i];
                int z = ztab[j];
                for (k = 0; k < isolength; k++)
                {
                    float newmz = mz + (isorange[k] / ((float)z));
                    // nearfast needs serious optimization or another way to do it... 
                    int pos = nearfast(dataMZ, newmz, lengthmz);
                    // reassignment might be able to be optimized. 
                    isotopepos[index3D(numz, isolength, i, j, k)] = pos;
                }
            }
        }
    }
    // definitely don't need two loops here. 
#pragma omp parallel for private (i,j,k), schedule(auto)
    for (i = 0; i < lengthmz; i++)
    {
        for (j = 0; j < numz; j++)
        {
            // again, if statement should be optimizable. 
            if (barr[index2D(numz, i, j)] == 1)
            {
                float mass = mtab[index2D(numz, i, j)];
                float mid = isotopemid(mass, isoparams);
                float sig = isotopesig(mass, isoparams);

                if (sig == 0) { printf("Error: Sigma Isotope Parameter is 0"); exit(102); }
                float alpha = isotopealpha(mass, isoparams);
                float amp = (1.0 - alpha) / (sig * 2.50662827);
                float beta = isotopebeta(mass, isoparams);
                float tot = 0;
                for (k = 0; k < isolength; k++)
                {
                    float e = alpha * exp(-isoindex[k] * beta);
                    float g = amp * exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
                    float temp = e + g;
                    tot += temp;
                    isotopeval[index3D(numz, isolength, i, j, k)] = temp;
                }
                for (k = 0; k < isolength; k++)
                {
                    if (tot > 0) { isotopeval[index3D(numz, isolength, i, j, k)] = isotopeval[index3D(numz, isolength, i, j, k)] / tot; }
                }
            }
        }
    }
    free(isorange);
    free(isoindex);
}

void isotope_dist(float mass, int isolength, int* isoindex, float* isovals, float* isoparams)
{
    float mid = isotopemid(mass, isoparams);
    float sig = isotopesig(mass, isoparams);
    if (sig == 0) { printf("Error: Sigma Isotope Parameter is 0"); exit(102); }
    float alpha = isotopealpha(mass, isoparams);
    float amp = 1.0 - alpha;
    float beta = isotopebeta(mass, isoparams);
    float tot = 0;
    int k;
    for (k = 0; k < isolength; k++)
    {
        float e = alpha * exp(-isoindex[k] * beta);
        float g = amp / (sig * 2.50662827) * exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
        tot += e + g;
        isovals[k] = e + g;
    }
    for (k = 0; k < isolength; k++)
    {
        if (tot > 0) { isovals[k] = isovals[k] / tot; }
    }
}


void test_isotopes(float mass, float* isoparams)
{
    int i;
    float maxmid = isotopemid(mass, isoparams);
    float maxsig = isotopesig(mass, isoparams);

    int isostart = 0;
    int isoend = (int)(maxmid + 4 * maxsig);
    if (isoend < 4) { isoend = 4; }
    int isolength = isoend - isostart;
    float* isorange = NULL;
    int* isoindex = NULL;
    isorange = calloc(isolength, sizeof(float));
    isoindex = calloc(isolength, sizeof(int));
    for (i = 0; i < isolength; i++)
    {
        isoindex[i] = (isostart + i);
    }
    isotope_dist(mass, isolength, isoindex, isorange, isoparams);
    free(isorange);
    free(isoindex);
}



void setup_and_make_isotopes(Config* config, Input* inp) {

    config->isolength = setup_isotopes(inp->isoparams, inp->isotopepos, inp->isotopeval, inp->mtab, inp->nztab, 
        inp->barr, inp->dataMZ, config->lengthmz, config->numz);

    inp->isotopepos = calloc(config->isolength * config->lengthmz * config->numz, sizeof(int));
    inp->isotopeval = calloc(config->isolength * config->lengthmz * config->numz, sizeof(float));

    make_isotopes(inp->isoparams, inp->isotopepos, inp->isotopeval, inp->mtab, inp->nztab, inp->barr, inp->dataMZ, config->lengthmz, config->numz);

}
void SetupAndMakeIsotopes(Config config, Input inp) {
    config.isolength = setup_isotopes(config.isolength, inp.isotopepos, inp.isotopeval, inp.mtab, inp.nztab,
        inp.barr, inp.dataMZ, config.lengthmz, config.numz); 
    inp.isotopepos = calloc(config.isolength * config.lengthmz * config.numz, sizeof(int)); 
    inp.isotopeval = calloc(config.isolength * config.lengthmz * config.numz, sizeof(float)); 

    make_isotopes(inp.isoparams, inp.isotopepos, inp.isotopeval, inp.mtab, inp.nztab, inp.barr, inp.dataMZ,
        config.lengthmz, config.numz); 
}
void monotopic_to_average(const int lengthmz, const int numz, float* blur, const char* barr, int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval)
{
    float* newblur = NULL;
    newblur = calloc(lengthmz * numz, sizeof(float));
    unsigned int i, j, k;
    for (i = 0; i < lengthmz; i++)
    {
        for (j = 0; j < numz; j++)
        {
            if (barr[index2D(numz, i, j)] == 1) {
                float topval = blur[index2D(numz, i, j)];
                for (k = 0; k < isolength; k++)
                {
                    int pos = isotopepos[index3D(numz, isolength, i, j, k)];
                    float val = isotopeval[index3D(numz, isolength, i, j, k)];
                    newblur[index2D(numz, pos, j)] += topval * (float)val;
                }
            }
        }
    }
    memcpy(blur, newblur, sizeof(float) * lengthmz * numz);
    free(newblur);
}

// Isotopes //




