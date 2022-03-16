//
//  UniDec_Main.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "UniDec_Main.h"
void ReadInputs(Input *inp, Config *config) {
    inp->nztab = calloc(config->numz, sizeof(int));
    inp->mtab = calloc(config->lengthmz * config->numz, sizeof(float)); 
    inp->barr = calloc(config->lengthmz * config->numz, sizeof(float)); 
    for (int i = 0; i < config->numz; i++) {
        inp->nztab[i] = i + config->startz;
    }
    for (int j = 0; j < config->numz; j++)
    {
        if (inp->nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
    }
    //Test to make sure no two data points has the same x value
    for (int i = 0; i < config->lengthmz - 1; i++)
    {
        if (inp->dataMZ[i] == inp->dataMZ[i + 1]) { printf("Error: Two data points are identical: %f %f \n\n", inp->dataMZ[i], inp->dataMZ[i + 1]); exit(104); }
    }
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < config->lengthmz; i++)
    {
        for (int j = 0; j < config->numz; j++)
        {
            inp->mtab[index2D(config->numz, i, j)] = (inp->dataMZ[i] * inp->nztab[j] - config->adductmass * inp->nztab[j]);
        }
    }
    ignorezeros(inp->barr, inp->dataInt, config->lengthmz, config->numz);
}
Input ReadInputsByValue(Input inp, Config* config) {
    inp.nztab = calloc(config->numz, sizeof(int));
    inp.mtab = calloc(config->lengthmz * config->numz, sizeof(float));
    inp.barr = calloc(config->lengthmz * config->numz, sizeof(float));
    for (int i = 0; i < config->numz; i++) {
        inp.nztab[i] = i + config->startz;
    }
    for (int j = 0; j < config->numz; j++)
    {
        if (inp.nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
    }
    //Test to make sure no two data points has the same x value
    for (int i = 0; i < config->lengthmz - 1; i++)
    {
        if (inp.dataMZ[i] == inp.dataMZ[i + 1]) { exit(104); }
    }
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < config->lengthmz; i++)
    {
        for (int j = 0; j < config->numz; j++)
        {
            inp.mtab[index2D(config->numz, i, j)] = (inp.dataMZ[i] * inp.nztab[j] - config->adductmass * inp.nztab[j]);
        }
    }
    ignorezeros(inp.barr, inp.dataInt, config->lengthmz, config->numz);
    return inp; 
}


void SetLimits(Config config, Input* inp) {
    //Determines the indexes of each test mass from mfile in m/z space
    int* testmasspos = malloc(sizeof(float) * config.mfilelen * config.numz);
    if (config.mflag == 1 && config.limitflag == 1) {
        for (int i = 0; i < config.mfilelen; i++)
        {
            for (int j = 0; j < config.numz; j++)
            {
                float mztest = (inp->testmasses[i] + config.adductmass * inp->nztab[j]) / (float)inp->nztab[j];
                testmasspos[index2D(config.numz, i, j)] = nearfast(inp->dataMZ, mztest, config.lengthmz);
            }
        }
    }

    //If there is a mass file read, it will only allow masses close to those masses within some config.mtabsig window.
    if (config.mflag == 1 && config.limitflag == 0)
    {
        TestMassListWindowed(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab, inp->testmasses, config.mfilelen, config.mtabsig);
    }
    //If there is a mass file read and the mass table window (config.mtabsig) is 0, it will only write intensities at the m/z values closest to the m/z values read in from the mfile.
    else if (config.mflag == 1 && config.limitflag == 1)
    {
        TestMassListLimit(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab, testmasspos, config.mfilelen);
    }
    //Normally, write the intensity values if the values fall within the mass upperbound and lower bound
    else
    {
        TestMass(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab);
    }
    free(testmasspos);
}


Decon run_unidec(Input inp, Config config) {
    
    inp = SetupInputs();
    ReadInputs(&inp, &config);
    //Sets limits based on mass range and any test masses
    SetLimits(config, &inp);

    //Setup Isotope Distributions
    if (config.isotopemode > 0)
    {
        setup_and_make_isotopes(&config, &inp);
    }

    //................................................................
    //
    // Deconvolution
    //
    //...................................................................

    //Setup the Deconvolution
    Decon decon = SetupDecon();
    //Run the main Deconvolution		
    decon = MainDeconvolution(config, inp);
    return decon; 
}


Decon MainDeconvolution(const Config config, const Input inp)
{
    // get rid of all the file import 
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

    //...................................................................
    //
    //     Sets the mzdist with the peak shape
    //
    //....................................................................
    int numberOfElementsInBarr = config.lengthmz * config.numz; 
    barr = calloc(numberOfElementsInBarr, sizeof(char));
    memcpy(barr, inp.barr, sizeof(barr));

    //Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
    float threshold = config.psthresh * fabs(config.mzsig) * config.peakshapeinflate;
    //Create a list of start and end values to box in arrays based on the above threshold
    starttab = calloc(config.lengthmz, sizeof(int));
    endtab = calloc(config.lengthmz, sizeof(int));
    int maxlength = 1;
    if (config.mzsig != 0) {
        //Gets maxlength and sets start and endtab
        maxlength = SetStartsEnds(config, &inp, starttab, endtab, threshold);

        //Changes dimensions of the peak shape function. 1D for speedy and 2D otherwise
        int pslen = config.lengthmz;
        if (config.speedyflag == 0) { pslen = config.lengthmz * maxlength; }
   
        mzdist = calloc(pslen, sizeof(float)); //malloc(sizeof(float) * pslen); //
        //if (verbose == 1) { printf("mzdist: %p %p\n", mzdist, NULL); }
        //memset(mzdist, 0, pslen * sizeof(float));
        
        if (pslen * sizeof(float) / 1000000000 > 4) { printf("Danger: Your data may crash the memory. Consider setting the Peak FWHM to 0.\n"); }
        int makereverse = 0;
        if (config.mzsig < 0 || config.beta < 0) { makereverse = 1; rmzdist = calloc(pslen, sizeof(float));
        //memset(rmzdist, 0, pslen * sizeof(float));
        }
        else { rmzdist = calloc(0, sizeof(float)); }

        //Calculates the distance between mz values as a 2D or 3D matrix

        if (config.speedyflag == 0)
        {
           
            MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab, inp.dataMZ, fabs(config.mzsig) * config.peakshapeinflate, config.psfun, config.speedyflag, mzdist, rmzdist, makereverse);
        }
        else
        {
            
            //Calculates peak shape as a 1D list centered at the first element for circular convolutions
            MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag, fabs(config.mzsig) * config.peakshapeinflate, config.psfun, mzdist, rmzdist, makereverse);
        }

    }
    else
    {
        mzdist = calloc(0, sizeof(float));
        maxlength = 0;
    }

    //....................................................
    //
    //    Setting up the neighborhood blur
    //
    //......................................................

    //sets some parameters regarding the neighborhood blur function
    if (config.zsig >= 0 && config.msig >= 0) {
        zlength = 1 + 2 * (int)config.zsig;
        mlength = 1 + 2 * (int)config.msig;
    }
    else {
        if (config.zsig != 0) { zlength = 1 + 2 * (int)(3 * fabs(config.zsig) + 0.5); }
        else { zlength = 1; }
        if (config.msig != 0) { mlength = 1 + 2 * (int)(3 * fabs(config.msig) + 0.5); }
        else { mlength = 1; }
    }
    numclose = mlength * zlength;

    //Sets up the blur function in oligomer mass and charge
    mind = calloc(mlength, sizeof(int));
    mdist = calloc(mlength, sizeof(float));

    for (int i = 0; i < mlength; i++)
    {
        mind[i] = i - (mlength - 1) / 2;
        if (config.msig != 0) { mdist[i] = exp(-(pow((i - (mlength - 1) / 2.), 2)) / (2.0 * config.msig * config.msig)); }
        else { mdist[i] = 1; }
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

    //Initializing memory
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

    //Set up blur
    //MakeBlur(config.lengthmz, config.numz, numclose, barr, closezind, closemind, inp.mtab, config.molig, config.adductmass, inp.nztab, inp.dataMZ, closeind, threshold, config);
    MakeSparseBlur(numclose, barr, closezind, closemind, inp.mtab, inp.nztab, inp.dataMZ, closeind, closeval, closearray, config);

    int badness = 1;
    for (int i = 0; i < config.lengthmz * config.numz; i++)
    {
        if (barr[i] == 1) { badness = 0;}
    }
    if (badness == 1) { printf("ERROR: Setup is bad. No points are allowed\n"); exit(10); }

    //Determine the maximum intensity in the data
    float dmax = Max(inp.dataInt, config.lengthmz);
    float betafactor = 1;
    if (dmax > 1) { betafactor=dmax; }

    //...................................................
    //
    //  Setting up and running the iteration
    //
    //.........................................................


    //Applies the intensity threshold to kill peaks
    if (config.intthresh != -1) { KillB(inp.dataInt, barr, config.intthresh, config.lengthmz, config.numz, config.isolength, inp.isotopepos, inp.isotopeval); }

    //Creates an intial probability matrix, decon.blur, of 1 for each element
    decon.blur = calloc(config.lengthmz * config.numz, sizeof(float));
    decon.newblur = calloc(config.lengthmz * config.numz, sizeof(float));
    oldblur = calloc(config.lengthmz * config.numz, sizeof(float));

    if (config.baselineflag == 1) {
        printf("Auto Baseline Mode On: %d\n", config.aggressiveflag);
        decon.baseline = calloc(config.lengthmz, sizeof(float));
        decon.noise = calloc(config.lengthmz, sizeof(float));
    }


    //#pragma omp parallel for private (i,j), schedule(auto)
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
    if (config.baselineflag == 1)
    {
        if (config.mzsig != 0)
        {
            deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabs(config.mzsig));
            if (config.aggressiveflag == 2)
            {
                for (int i = 0; i < 10; i++)
                {
                    deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabs(config.mzsig));
                }
                for (int i = 0; i < config.lengthmz; i++)
                {
                    if (decon.baseline[i] > 0) {
                        dataInt2[i] -= decon.baseline[i];
                    }
                }
                //memcpy(decon.baseline, dataInt2, sizeof(float)*config.lengthmz);
            }
        }
        else
        {
            printf("Ignoring baseline subtraction because peak width is 0\n");
        }

    }


    //Run the iteration
    float blurmax = 0;
    decon.conv = 0;
    int off = 0;

    for (int iterations = 0; iterations < abs(config.numit); iterations++)
    {
        decon.iterations = iterations;
        if (config.beta > 0 && iterations > 0)
        {

            softargmax(decon.blur, config.lengthmz, config.numz, config.beta/betafactor);
            //printf("Beta %f\n", beta);
        }
        else if (config.beta < 0 && iterations >0)
        {
            softargmax_transposed(decon.blur, config.lengthmz, config.numz, fabs(config.beta/betafactor), barr, maxlength, config.isolength, inp.isotopepos, inp.isotopeval, config.speedyflag, starttab, endtab, rmzdist, config.mzsig);
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
                    diff += pow((decon.blur[i] - oldblur[i]), 2);
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
    decon.error = errfunspeedy(config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, &decon.rsquared);

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
                SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
                    inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
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
    decon.uniscore = score(config, &decon, inp, scorethreshold);

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

// got rid of autotune and run unidec. Focusing only on the main deconvolution function for mzLib import

