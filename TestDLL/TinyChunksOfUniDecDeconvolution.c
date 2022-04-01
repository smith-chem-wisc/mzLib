// Goal is to create tiny little chunks of UniDec's algorithm so that each component is completely tested and functional

#include "TinyChunksOfUniDecDeconvolution.h"

int MemoryAllocationOfBarr(void) {
    char* barr = NULL;

    int mlength, zlength, numclose,
        * mind = NULL,
        * zind = NULL,
        * closemind = NULL,
        * closezind = NULL,
        * closeind = NULL,
        * starttab = NULL,
        * endtab = NULL;
    float
        * mdist = NULL,
        * zdist = NULL,
        * mzdist = NULL,
        * rmzdist = NULL,
        * oldblur = NULL,
        * closeval = NULL,
        * closearray = NULL;

    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    free(barr);
    free(closeind);
    return 1; 
}

int AllocateMemoryToPointersThenFree() {
    char* barr = NULL;

    int mlength, zlength, numclose,
        * mind = calloc(10 * 10, sizeof(int)); 
        int * zind = NULL,
        * closemind = NULL,
        * closezind = NULL,
        * closeind = NULL,
        * starttab = NULL,
        * endtab = NULL;
    float
        * mdist = NULL,
        * zdist = NULL,
        * mzdist = NULL,
        * rmzdist = NULL,
        * oldblur = NULL,
        * closeval = NULL,
        * closearray = NULL;

    if (mind != NULL) {

        for (int i = 0; i < 10; i++) {
            mind[i] = i;
        }
    }

    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    free(barr);
    free(closeind);
    return 1; 
}

int UseConfigAndInputToCreatePointerValues(Input inp, Config config) {

}

char UseMemcpyInC() {
    char testBarr[2] = "12"; 
    char* barr = calloc(4, sizeof(char)); 
    memcpy(barr, testBarr, 4 * sizeof(char)); 
    return barr[0]; 
}

char UseMemcpyWithInpAndConfigObjects() {
    Config config; 
    Input inp; 

    float xarray[3] = { (float)1, (float)2, (float)3 }; 
    float yarray[3] = { (float)10, (float)20, (float)30 }; 
    char xBarr[3] = "001"; 

    inp.dataMZ = xarray; 
    inp.dataInt = yarray; 
    inp.barr = xBarr; 

    char* barr = calloc(3, sizeof(char)); 
    memcpy(barr, inp.barr, 3 * sizeof(char)); 


    return barr[2]; 
    
}

int MemoryObjectAllocationToHeap(Config config, Input inp) {
    // set up the memroy
    // pointers are set to null because not all of them are used. 
    Decon decon = SetupDecon(); 
    int numberOfBarrElements = 100; 

    //char* barr = calloc(numberOfBarrElements, sizeof(char)); 
    int* starttab = calloc(config.lengthmz, sizeof(int)); 
    int* endtab = calloc(config.lengthmz, sizeof(int));
    int zlength = 1 + 2 * (int)config.zsig; 
    int mlength = 1 + 2 * (int)config.msig; 
    int pslen = config.lengthmz; 
    float* mzdist = calloc(pslen, sizeof(float));
    int* mind = calloc(mlength, sizeof(int));
    
    float* mdist = calloc(mlength, sizeof(float));
    
    int* zind = calloc(zlength, sizeof(int));
    
    
    float* zdist = calloc(zlength, sizeof(float));
    int numclose = mlength * zlength; 
    int* closemind = calloc(numclose, sizeof(int));
    
    int* closezind = calloc(numclose, sizeof(int));
    int* closeval = calloc(numclose, sizeof(float));
    int* closeind = calloc(numclose * config.lengthmz * config.numz, sizeof(int));
    int* closearray = calloc(numclose * config.lengthmz * config.numz, sizeof(float));
    float* rmzdist = calloc(pslen, sizeof(float)); 
    decon.blur = calloc(config.lengthmz * config.numz, sizeof(float));
    decon.newblur = calloc(config.lengthmz * config.numz, sizeof(float));
    float* oldblur = calloc(config.lengthmz * config.numz, sizeof(float));
    decon.baseline = calloc(config.lengthmz, sizeof(float));
    decon.noise = calloc(config.lengthmz, sizeof(float));
    float* dataInt2 = calloc(config.lengthmz, sizeof(float));
    decon.fitdat = calloc(config.lengthmz, sizeof(float));
    
    float threshold = config.psthresh * (float)fabs((double)config.mzsig) * config.peakshapeinflate; 
    int maxlength = SetStartsEnds(config, &inp, starttab, endtab, threshold);
    //decon.error = errfunspeedy(config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, &decon.rsquared);

    memcpy(decon.newblur, decon.blur, (config.lengthmz * config.numz) * sizeof(float));
    decon.massaxis = calloc(decon.mlen, sizeof(float));
    decon.massaxisval = calloc(decon.mlen, sizeof(float));
    decon.massgrid = calloc(decon.mlen * config.numz, sizeof(float));
    
    float massmax = config.masslb; 
    float massmin = config.massub; 
    decon.mlen = (int)((massmax - massmin) / config.massbins); 
    //memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
    //memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(float));

    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    //free(barr);
    free(closeind);
    return 1; 
}

/*int MemoryObjectAllocationToHeapConfigPtr(Config* config, Input inp) {
    // Memory allocation to decon is the issue. 
    // set up the memroy
    // pointers are set to null because not all of them are used. 
    Decon decon = SetupDecon();
    int numberOfElementsInBarr = config->lengthmz * config->numz; 
    char* barr = calloc(numberOfElementsInBarr, sizeof(char)); 
    memcpy(barr, inp.barr, numberOfElementsInBarr * sizeof(char)); 
    int* starttab = calloc(config->lengthmz, sizeof(int));
    int* endtab = calloc(config->lengthmz, sizeof(int));
    int zlength = 1 + 2 * (int)config->zsig;
    int mlength = 1 + 2 * (int)config->msig;
    int pslen = config->lengthmz;
    float* mzdist = calloc(pslen, sizeof(float));
    int* mind = calloc(mlength, sizeof(int));

    float* mdist = calloc(mlength, sizeof(float));

    int* zind = calloc(zlength, sizeof(int));


    float* zdist = calloc(zlength, sizeof(float));
    int numclose = mlength * zlength;
    int* closemind = calloc(numclose, sizeof(int));

    int* closezind = calloc(numclose, sizeof(int));
    int* closeval = calloc(numclose, sizeof(float));
    int* closeind = calloc(numclose * config->lengthmz * config->numz, sizeof(int));
    int* closearray = calloc(numclose * config->lengthmz * config->numz, sizeof(float));
    float* rmzdist = calloc(pslen, sizeof(float));
    //decon.blur = calloc(config->lengthmz * config->numz, sizeof(float));
    //decon.newblur = calloc(config->lengthmz * config->numz, sizeof(float));
    float* oldblur = calloc(config->lengthmz * config->numz, sizeof(float));
    //decon.baseline = calloc(config->lengthmz, sizeof(float));
    //decon.noise = calloc(config->lengthmz, sizeof(float));
    float* dataInt2 = calloc(config->lengthmz, sizeof(float));
    //decon.fitdat = calloc(config->lengthmz, sizeof(float));

    float threshold = config->psthresh * (float)fabs((double)config->mzsig) * config->peakshapeinflate;
    int maxlength = SetStartsEnds(*config, &inp, starttab, endtab, threshold);
    decon.error = errfunspeedy(*config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, decon.rsquared);
    decon.error = errfunspeedy(*config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, decon.rsquared);

    //memcpy(decon.newblur, decon.blur, (config->lengthmz * config->numz) * sizeof(float));
    //decon.massaxis = calloc(decon.mlen, sizeof(float));
    //decon.massaxisval = calloc(decon.mlen, sizeof(float));
    decon.massgrid = calloc(decon.mlen * config->numz, sizeof(float));

    float massmax = config->masslb;
    float massmin = config->massub;
    //decon.mlen = (int)((massmax - massmin) / config->massbins);
    //memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
    memset(decon.massgrid, 0, decon.mlen * config->numz * sizeof(float));

    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    free(barr);
    free(closeind);
    FreeDecon(decon); 
    return 1;
}
*/
int TestSetStartEnds(Input inp, Config config) {
    int* starttab = calloc(config.lengthmz, sizeof(int));
    int* endtab = calloc(config.lengthmz, sizeof(int));
    float threshold = config.psthresh * fabs(config.mzsig) * config.peakshapeinflate;
    return SetStartsEnds(config, &inp, starttab, endtab, threshold); 
}
int TestFreeDecon() {
    Decon decon = SetupDecon(); 
    FreeDecon(decon); 
    return 1;
}
int TestSetupAndAllocateMemoryToDecon() {
    Decon decon = SetupDecon(); 
    decon.fitdat = calloc(100, sizeof(float)); 
    memset(decon.massaxisval, (float)0, decon.mlen * 100 * (int)sizeof(float)); 
    decon.blur = calloc(100*200, sizeof(float));
    decon.newblur = calloc(100*200, sizeof(float));
    decon.baseline = calloc(100*200, sizeof(float));
    decon.noise = calloc(100, sizeof(float));
    decon.fitdat = calloc(100, sizeof(float));
    // decon.error = errfunspeedy(*config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, &decon.rsquared);
    memcpy(decon.newblur, decon.blur, 100 * sizeof(float));
    decon.massaxis = calloc(decon.mlen, sizeof(float));
    decon.massaxisval = calloc(decon.mlen, sizeof(float));
    decon.massgrid = calloc(decon.mlen * 100, sizeof(float));
    decon.mlen = (int)((100) / 10);
    memset(decon.massaxisval, 0, 100 * sizeof(float)); 
    memset(decon.massgrid, 0, 100 * sizeof(float)); 
    FreeDecon(decon);
    return 1; 
}


Decon TestSetupAndReturnDecon() {
    Decon decon = SetupDecon();    
    decon.massaxis = calloc(10, sizeof(float));
    return decon;
}

Decon MainDeconWithMinimalControlFlow(Config config, Input inp) {
    Decon decon = SetupDecon();
    char* barr = NULL;

    int mlength, zlength, numclose,
        * mind = NULL,
        * zind = NULL,
        * closemind = NULL,
        * closezind = NULL,
        * closeind = NULL,
        * starttab = NULL,
        * endtab = NULL;
    float
        * mdist = NULL,
        * zdist = NULL,
        * mzdist = NULL,
        * rmzdist = NULL,
        * oldblur = NULL,
        * closeval = NULL,
        * closearray = NULL;

    int numberOfElementsInBarr = config.lengthmz * config.numz;
    barr = calloc(numberOfElementsInBarr, sizeof(char));
    memcpy(barr, inp.barr, sizeof(barr));
    //Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
    float threshold = config.psthresh * fabs(config.mzsig) * config.peakshapeinflate;
    //Create a list of start and end values to box in arrays based on the above threshold
    starttab = calloc(config.lengthmz, sizeof(int));
    endtab = calloc(config.lengthmz, sizeof(int));
    int maxlength = 1;
    maxlength = SetStartsEnds(config, &inp, starttab, endtab, threshold);
    int pslen = config.lengthmz * maxlength;
    mzdist = calloc(pslen, sizeof(float));
    

    //Sets up the blur function in oligomer mass and charge


    int makereverse = 1; 
    rmzdist = calloc(pslen, sizeof(float));
    // fabs takes returns the abosolute values of that number. 
    MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab, inp.dataMZ, fabs(config.mzsig) * config.peakshapeinflate, config.psfun, config.speedyflag, mzdist, rmzdist, makereverse);
    zlength = 1 + 2 * (int)config.zsig;
    mlength = 1 + 2 * (int)config.msig;
    mind = calloc(mlength, sizeof(int));
    mdist = calloc(mlength, sizeof(float));
    
    if (mind != NULL) {
        for (int i = 0; i < mlength; i++)
        {
            mind[i] = i - (mlength - 1) / 2;
            if (config.msig != 0) { mdist[i] = exp(-(pow((i - (mlength - 1) / 2.), 2)) / (2.0 * config.msig * config.msig)); }
            else { mdist[i] = 1; }
        }
    }
    
    
    zind = calloc(zlength, sizeof(int));
    zdist = calloc(zlength, sizeof(float));
    for (int i = 0; i < zlength; i++)
    {
        zind[i] = i - (zlength - 1) / 2;
        if (config.zsig != 0) { zdist[i] = exp(-(pow((i - (zlength - 1) / 2.), 2)) / (2.0 * config.zsig * config.zsig)); }
        else { zdist[i] = 1; }
        //printf("%f\n", zdist[i]);
    }

    numclose = mlength * zlength;
    closemind = calloc(numclose, sizeof(int));
    closezind = calloc(numclose, sizeof(int)); 
    closeval = calloc(numclose, sizeof(float));
    closeind = calloc(numclose * config.lengthmz * config.numz, sizeof(int));
    closearray = calloc(numclose * config.lengthmz * config.numz, sizeof(float));


    //Determines the indexes of things that are close as well as the values used in the neighborhood convolution
    for (int k = 0; k < numclose; k++)
    {
        closemind[k] = mind[k % mlength];
        closezind[k] = zind[(int)k / mlength];
        closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
    }
    simp_norm_sum(mlength, mdist);
    simp_norm_sum(zlength, zdist);
    simp_norm_sum(numclose, closeval);

    MakeSparseBlur(numclose, barr, closezind, closemind, inp.mtab, inp.nztab, inp.dataMZ, closeind, closeval, closearray, config);

    int badness = 1;
    for (int i = 0; i < config.lengthmz * config.numz; i++)
    {
        if (barr[i] == 1) { badness = 0; }
    }
    if (badness == 1) { printf("ERROR: Setup is bad. No points are allowed\n"); exit(10); }

    //Determine the maximum intensity in the data
    float dmax = Max(inp.dataInt, config.lengthmz);
    float betafactor = 1;
    if (dmax > 1) { betafactor = dmax; }
    
    //Something wrong with KillB. 
    TestingKillBFunction(inp.dataInt, barr, config.intthresh, config.lengthmz, config.numz, 
        config.isolength, inp.isotopepos, inp.isotopeval);
    decon.blur = calloc(config.lengthmz * config.numz, sizeof(float));
    decon.newblur = calloc(config.lengthmz * config.numz, sizeof(float));
    oldblur = calloc(config.lengthmz * config.numz, sizeof(float));

    decon.baseline = calloc(config.lengthmz, sizeof(float));
    decon.noise = calloc(config.lengthmz, sizeof(float));

    // Create decon blur
    for (int i = 0; i < config.lengthmz; i++)
    {
        float val = inp.dataInt[i] / ((float)(config.numz + 2));
        if (config.baselineflag == 1) {

            decon.baseline[i] = val;
            decon.noise[i] = val;
        }

        for (int j = 0; j < config.numz; j++)
        {
            if (barr[index2D(config.numz, i, j)] == 1) {
                if (config.isotopemode == 0) {
                    decon.blur[index2D(config.numz, i, j)] = val;
                }
                else { decon.blur[index2D(config.numz, i, j)] = 1; }
            }
            else
            {
                decon.blur[index2D(config.numz, i, j)] = 0;
            }
        }
    }

    memcpy(oldblur, decon.blur, sizeof(float) * config.lengthmz * config.numz);
    memcpy(decon.newblur, decon.blur, sizeof(float) * config.lengthmz * config.numz);

    float* dataInt2 = NULL;
    dataInt2 = calloc(config.lengthmz, sizeof(float));
    memcpy(dataInt2, inp.dataInt, sizeof(float) * config.lengthmz);
    deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabs(config.mzsig));
    float blurmax = 0;
    decon.conv = 0;
    int off = 0;

    for (int iterations = 0; iterations < abs(config.numit); iterations++)
    {
        decon.iterations = iterations;
        if (config.beta > 0 && iterations > 0)
        {

            softargmax(decon.blur, config.lengthmz, config.numz, config.beta / betafactor);
            //printf("Beta %f\n", beta);
        }
        else if (config.beta < 0 && iterations >0)
        {
            softargmax_transposed(decon.blur, config.lengthmz, config.numz, fabs(config.beta / betafactor), barr, maxlength, config.isolength, inp.isotopepos, inp.isotopeval, config.speedyflag, starttab, endtab, rmzdist, config.mzsig);
        }

        if (config.psig >= 1 && iterations > 0)
        {
            point_smoothing(decon.blur, barr, config.lengthmz, config.numz, abs((int)config.psig));
            //printf("Point Smoothed %f\n", config.psig);
        }
        else if (config.psig < 0 && iterations >0)
        {
            point_smoothing_peak_width(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, config.speedyflag, barr);
        }


        //Run Blurs
        if (config.zsig >= 0 && config.msig >= 0) {
            blur_it_mean(config.lengthmz, config.numz, numclose, closeind, decon.newblur, decon.blur, barr, closearray, config.zerolog);
        }
        else if (config.zsig > 0 && config.msig < 0)
        {
            blur_it_hybrid1(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
        }
        else if (config.zsig < 0 && config.msig > 0)
        {
            blur_it_hybrid2(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
        }
        else {
            blur_it(config.lengthmz, config.numz, numclose, closeind, closearray, decon.newblur, decon.blur, barr);
        }

        //Run Richardson-Lucy Deconvolution
        deconvolve_iteration_speedy(config.lengthmz, config.numz, maxlength,
            decon.newblur, decon.blur, barr, config.aggressiveflag, dataInt2,
            config.isolength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, rmzdist, config.speedyflag,
            config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);

        //Determine the metrics for conversion. Only do this every 10% to speed up.
        if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations>0.9 * config.numit)) {
            float diff = 0;
            float tot = 0;
            for (int i = 0; i < config.lengthmz * config.numz; i++)
            {
                if (barr[i] == 1)
                {
                    diff += pow(((double)decon.blur[i] - (double)oldblur[i]), 2);
                    tot += decon.blur[i];
                }
            }
            if (tot != 0) { decon.conv = (diff / tot); }
            else { decon.conv = 12345678; printf("m/z vs. charge grid is zero. Iteration: %d\n", iterations); }

            //printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
            if (decon.conv < 0.000001) {
                if (off == 1 && config.numit > 0) {
                    printf("Converged in %d iterations.\n\n", iterations);
                    break;
                }
                off = 1;
            }
            memcpy(oldblur, decon.blur, config.lengthmz * config.numz * sizeof(float));
        }

    }
    free(dataInt2);
    free(oldblur);

    //................................................................
    //
    //     Setting up the outputs
    //
    //...............................................................



    //Reset the peak shape if it was inflated
    if (config.peakshapeinflate != 1 && config.mzsig != 0) {
        if (config.speedyflag == 0)
        {
            MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab, inp.dataMZ, fabs(config.mzsig), config.psfun, config.speedyflag, mzdist, rmzdist, 0);
        }
        else
        {
            MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag, fabs(config.mzsig), config.psfun, mzdist, rmzdist, 0);
        }
        printf("mzdist reset: %f\n", config.mzsig);
    }


    //Determine the maximum intensity in the blur matrix
    blurmax = Max(decon.blur, config.lengthmz * config.numz);
    float cutoff = 0;
    if (blurmax != 0) { cutoff = 0.000001; }

    //Apply The Cutoff
    ApplyCutoff1D(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);


    //Calculate the fit data and error.
    decon.fitdat = calloc(config.lengthmz, sizeof(float));
    /*decon.error = errfunspeedy(config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval,
        starttab, endtab, mzdist, &decon.rsquared);*/

    //Fix issues with fitdat and consecutive zero data points
    //TODO: It might be possible to build this in to convolve_simp so that this isn't necessary but it would require a 1D barr.
    if (config.intthresh != -1)
    {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < config.lengthmz-1; i++)
        {
            if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
            {
                decon.fitdat[i] = 0;
                decon.fitdat[i + 1] = 0;
            }
        }
    }


    // Charge scaling (orbimode)
    if (config.orbimode == 1)
    {
        printf("Rescaling charge states and normalizing ");
        charge_scaling(decon.blur, inp.nztab, config.lengthmz, config.numz);
        //simp_norm(config.lengthmz * config.numz, decon.blur);
        printf("Done\n");
    }

    //Change Monoisotopic to Average if necessary
    if (config.isotopemode == 2)
    {
        monotopic_to_average(config.lengthmz, config.numz, decon.blur, barr, config.isolength, inp.isotopepos, inp.isotopeval);
    }

    //newblur is repurposed as the convolution of blur by the mz peaks shape
    float newblurmax = blurmax;
    if ((config.rawflag == 0 || config.rawflag == 2)) {
        if (config.mzsig != 0) {
            newblurmax = Reconvolve(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, decon.newblur, config.speedyflag, barr);
        }
        else
        {
            memcpy(decon.newblur, decon.blur, (config.lengthmz * config.numz) * sizeof(float));
        }
    }


    //.......................................................
    //
    //  Mass space outputs
    //
    //..........................................................

    //Determine the maximum and minimum allowed masses.
    float massmax = config.masslb;
    float massmin = config.massub;
    if (config.fixedmassaxis == 0) {
        for (int i = 0; i < config.lengthmz; i++)
        {
            for (int j = 0; j < config.numz; j++)
            {
                if (decon.newblur[index2D(config.numz, i, j)] * barr[index2D(config.numz, i, j)] > newblurmax * cutoff)
                {
                    float testmax = inp.mtab[index2D(config.numz, i, j)] + threshold * inp.nztab[j]+config.massbins;
                    float testmin = inp.mtab[index2D(config.numz, i, j)] - threshold * inp.nztab[j];

                    //To prevent really wierd decimals
                    testmin = round(testmin / config.massbins) * config.massbins;
                    testmax = round(testmax / config.massbins) * config.massbins;

                    if (testmax > massmax){    massmax = testmax;}
                    if (testmin < massmin){    massmin = testmin;}
                }
            }
        }
        printf("Massmin: %f  ", massmin);
        printf("Massmax: %f  ", massmax);
    }
    else { massmax = config.massub; massmin = config.masslb; }

    //Checks to make sure the mass axis is good and makes a dummy axis if not
    decon.mlen = (int)(massmax - massmin) / config.massbins;
    if (decon.mlen < 1) {
        printf("Bad mass axis length: %d\n", decon.mlen);
        massmax = config.massub;
        massmin = config.masslb;
        decon.mlen = (int)(massmax - massmin) / config.massbins;

        //Declare the memory
        decon.massaxis = calloc(decon.mlen, sizeof(float));
        decon.massaxisval = calloc(decon.mlen, sizeof(float));
        decon.massgrid = calloc(decon.mlen * config.numz, sizeof(float));
        memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
        memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(float));

        //Create the mass axis
        for (int i = 0; i < decon.mlen; i++)
        {
            decon.massaxis[i] = massmin + i * config.massbins;
        }
        decon.uniscore = 0;
        printf("ERROR: No masses detected.\n");
    }
    else {

        //Declare the memory
        decon.massaxis = calloc(decon.mlen, sizeof(float));
        decon.massaxisval = calloc(decon.mlen, sizeof(float));
        decon.massgrid = calloc(decon.mlen * config.numz, sizeof(float));
        memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
        memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(float));


        //Create the mass axis
        for (int i = 0; i < decon.mlen; i++)
        {
            decon.massaxis[i] = massmin + i * config.massbins;
        }

        //Determine the mass intensities from m/z grid
        if (config.poolflag == 0) {
            if (config.rawflag == 1 || config.rawflag == 3) {
                IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.blur, decon.massgrid);
            }
            if (config.rawflag == 0 || config.rawflag == 2) {
                IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.newblur, decon.massgrid);
            }
        }
        else if (config.poolflag == 1) {
            if (config.rawflag == 1 || config.rawflag == 3) {
                InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
                    inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
            }
            if (config.rawflag == 0 || config.rawflag == 2) {
                InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
                    inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
            }
        }
        else if (config.poolflag == 2) {
            if (config.rawflag == 1 || config.rawflag == 3) {
                SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, 
                    decon.massaxis, config.adductmass, inp.dataMZ, decon.massgrid, 
                    decon.massaxisval, decon.blur);
            }
            if (config.rawflag == 0 || config.rawflag == 2) {
                SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
                    inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
            }
        }
        else {
            printf("Invalid poolflag %d\n", config.poolflag);
            exit(1987);
        }

    //.....................................
    // Scores
    // .......................................

    //Note this will not execute if the mass axis is bad
    float scorethreshold = 0;
    decon.uniscore = performScoring(config, &decon, inp, scorethreshold);

    }
    //Free Memory
    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    free(barr);
    free(closeind);

    return decon; 
}

Decon RunUniDecWithTestMainDeconAlgo(Input inp, Config config) {
    // Called by C# calling code instead: 
    // ReadInputsByValue(inp, &config);
    //Sets limits based on mass range and any test masses

    SetLimits(config, &inp);
    
    //Setup Isotope Distributions
    if (config.isotopemode > 0)
    {
        //setup_and_make_isotopes(&config, &inp);
    }
    /*
    //................................................................
    //
    // Deconvolution
    //
    //...................................................................

    //Run the main Deconvolution		
    */
    Decon result = MainDeconWithMinimalControlFlow(config, inp);
    return result; 
}

void TestingKillBFunction(float* I, char* B, float intthresh, int lengthmz, int numz, const int isolength, int* isotopepos, float* isotopeval)
{
    unsigned int i, j, k;
    if (isolength == 0) {
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                if (I[i] <= intthresh) { B[index2D(numz, i, j)] = 0; }
            }
        }
    }
    else
    {
        float cutoff = 0.5;
        for (i = 0; i < lengthmz; i++)
        {
            for (j = 0; j < numz; j++)
            {
                float max = 0;
                for (k = 0; k < isolength; k++)
                {
                    float val = isotopeval[index3D(numz, isolength, i, j, k)];
                    if (val > max) { max = val; }
                    if (val > cutoff * max) {
                        int pos = isotopepos[index3D(numz, isolength, i, j, k)];
                        if (I[pos] <= intthresh) { B[index2D(numz, i, j)] = 0; }
                    }
                }

            }
        }
    }
}
void TestingCharArrayMarshalling(char* arrayOfChar, int lengthArray) {
    // modify char array.
    for (int i = 0; i < lengthArray; i++) {
        arrayOfChar[i] = '1'; 
    }
}
