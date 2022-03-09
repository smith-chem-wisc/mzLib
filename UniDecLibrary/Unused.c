//
//  Unused.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Unused.h"

//Print an int array
void IntPrint(const int* array, const int length)
{
    for (int i = 0; i < length; i++)
    {
        printf("%d\n", array[i]);
    }
}

//Print a float array
void floatPrint(const float* array, const int length)
{
    for (int i = 0; i < length; i++)
    {
        printf("%f\n", array[i]);
    }
}


void textvectorprint(float* arr, int length)
{
    int levels = 20;
    float grad = 1.0 / ((float)levels);
    int i, j;
    printf("\n");
    float max = 0;
    for (i = 0; i < length; i++) {
        if (arr[i] > max) { max = arr[i]; }
    }
    if (max != 0) {
        for (i = 0; i < length; i++) {
            arr[i] = arr[i] / max;
        }
    }
    for (i = 0; i < levels; i++) {
        for (j = 0; j < length; j++)
        {
            if (arr[j] > grad * (levels - i)) { printf("| "); }
            else { printf("  "); }
        }
        printf("\n");
    }
}


// Unused //

// "Returns" discrete Fourier transform of input array. Unused.
// Use fftw library instead
void discretefouriertransform(double* input, double** output, int length) {
    double pi = 3.1415926535;
    for (int i = 0; i < length; i++) {
        output[i][0] = 0.0; // Real term
        output[i][1] = 0.0; // Imaginary term
        for (int j = 0; j < length; j++) {
            double inner = 2.0 * pi / length * i * j;
            output[i][0] += input[j] * cos(inner);
            output[i][1] += input[j] * (-1.0) * sin(inner);
        }
    }
}

// "Returns" discrete inverse Fourier transform of input. Unused.
// Use fftw library instead
void inversefouriertransform(double** input, double* output, int length) {
    double pi = 3.1415926535;
    for (int i = 0; i < length; i++) {
        output[i] = 0.0;    // Real term
        double imag = 0.0;    // Imaginary term
        for (int j = 0; j < length; j++) {    // (ac - bd) + i(ad + bc)
            double inner = 2.0 * pi / length * i * j;
            double c = cos(inner);
            double d = sin(inner);
            output[i] += (input[j][0] * c) - (input[j][1] * d);
            imag += (input[j][0] * d) + (input[j][1] * c);
        }
        output[i] = output[i] / ((double)length);
        imag = imag / ((double)length);
        printf("Imaginary term is %.2f", imag);
    }
}


void WriteDecon(const Config config, const Decon* decon, const Input* inp)
{
    hid_t file_id = config.file_id;

    char outdat[1024];
    FILE* out_ptr = NULL;
    //Write the fit data to a file.
    if (config.rawflag >= 0) {
        if (config.filetype == 0) {
            char outstring2[500];
            sprintf(outstring2, "%s_fitdat.bin", config.outfile);
            out_ptr = fopen(outstring2, "wb");
            fwrite(decon->fitdat, sizeof(float), config.lengthmz, out_ptr);
            fclose(out_ptr);
            //printf("Fit: %s\t", outstring2);
        }
        else {
            //strjoin(dataset, "/fit_data", outdat);
            //mh5writefile1d(file_id, outdat, config.lengthmz, decon->fitdat);
        }
    }

    //Write the baseline to a file.
    if (config.rawflag >= 0 && config.baselineflag == 1) {
        if (config.filetype == 0) {
            char outstring2[500];
            sprintf(outstring2, "%s_baseline.bin", config.outfile);
            out_ptr = fopen(outstring2, "wb");
            fwrite(decon->baseline, sizeof(float), config.lengthmz, out_ptr);
            fclose(out_ptr);
            //printf("Background: %s\t", outstring2);
        }
        else {
            //strjoin(dataset, "/baseline", outdat);
            //mh5writefile1d(file_id, outdat, config.lengthmz, decon->baseline);//
        }
    }

    //Writes the convolved m/z grid in binary format
    if (config.rawflag == 0 || config.rawflag == 1)
    {
        // Note to self
        // rawflag=0 -> Reconvolved/Profile -> newblur
        // rawflag=1 -> Raw/Centroid -> blur

        if (config.filetype == 0) {
            char outstring9[500];
            sprintf(outstring9, "%s_grid.bin", config.outfile);
            out_ptr = fopen(outstring9, "wb");
            if (config.rawflag == 0) { fwrite(decon->newblur, sizeof(float), config.lengthmz * config.numz, out_ptr); }
            if (config.rawflag == 1) { fwrite(decon->blur, sizeof(float), config.lengthmz * config.numz, out_ptr); }
            fclose(out_ptr);
            //printf("m/z grid: %s\t", outstring9);
        }
        else {
            strjoin(config.dataset, "/mz_grid", outdat);
            if (config.rawflag == 0) { mh5writefile1d(file_id, outdat, config.lengthmz * config.numz, decon->newblur); }
            if (config.rawflag == 1) { mh5writefile1d(file_id, outdat, config.lengthmz * config.numz, decon->blur); }

            float* chargedat = NULL;
            chargedat = calloc(config.numz, sizeof(float));
            float* chargeaxis = NULL;
            chargeaxis = calloc(config.numz, sizeof(float));

            for (int j = 0; j < config.numz; j++) {
                float val = 0;
                chargeaxis[j] = (float)inp->nztab[j];
                for (int i = 0; i < config.lengthmz; i++) {

                    val += decon->newblur[index2D(config.numz, i, j)];
                }
                chargedat[j] = val;
            }
            strjoin(config.dataset, "/charge_data", outdat);
            mh5writefile2d(file_id, outdat, config.numz, chargeaxis, chargedat);
            free(chargedat);
            free(chargeaxis);
        }
    }
    else if (config.filetype == 1) {
        // If the rawfile flag not 1 or 2 and this is an HDF5 file, zero out these arrays.
        strjoin(config.dataset, "/mz_grid", outdat);
        delete_group(file_id, outdat);
        strjoin(config.dataset, "/charge_data", outdat);
        delete_group(file_id, outdat);
        strjoin(config.dataset, "/mass_grid", outdat);
        delete_group(file_id, outdat);
    }

    //Writes the convolved mass grid in binary format
    if (config.rawflag == 0 || config.rawflag == 1) {
        if (config.filetype == 0) {
            char outstring10[500];
            sprintf(outstring10, "%s_massgrid.bin", config.outfile);
            out_ptr = fopen(outstring10, "wb");
            fwrite(decon->massgrid, sizeof(float), decon->mlen * config.numz, out_ptr);
            fclose(out_ptr);
            //printf("Mass Grid: %s\t", outstring10);
        }
        else {
            strjoin(config.dataset, "/mass_grid", outdat);
            mh5writefile1d(file_id, outdat, decon->mlen * config.numz, decon->massgrid);
        }
    }

    //Writes the mass values convolved with the peak shape
    if (config.rawflag == 0 || config.rawflag == 1 || config.rawflag == 2 || config.rawflag == 3) {
        if (config.filetype == 0) {
            char outstring4[500];
            sprintf(outstring4, "%s_mass.txt", config.outfile);
            out_ptr = fopen(outstring4, "w");
            for (int i = 0; i < decon->mlen; i++)
            {
                fprintf(out_ptr, "%f %f\n", decon->massaxis[i], decon->massaxisval[i]);
            }
            fclose(out_ptr);
            //printf("Masses: %s\n", outstring4);
        }
        else {
            //int mlen = remove_middle_zeros(decon->massaxis, decon->massaxisval, decon->mlen); //Doesn't really work
            strjoin(config.dataset, "/mass_data", outdat);
            mh5writefile2d(file_id, outdat, decon->mlen, decon->massaxis, decon->massaxisval);
        }
    }

    if (config.filetype == 1 && decon->plen > 0) {
        WritePeaks(config, decon);
    }

}
// This function reads the H5 file and inputs the contents into the config and input structs. This function will be unecessary if directly passing from C# to the .dll.
// Need to refactor this code to convert mass spectra data to the Input object to pass directly to C. 
void ReadInputs(int argc, char* argv[], Config* config, Input* inp)
{
    if (config->filetype == 1) {
        if (config->metamode != -2)
        {
            strcpy(config->dataset, "/ms_dataset");
            char strval[1024];
            sprintf(strval, "/%d", config->metamode);
            strcat(config->dataset, strval);
            printf("HDF5 Data Set: %s\n", config->dataset);
        }
        else
        {
            strcpy(config->dataset, "/ms_data");
        }

        char outdat[1024];
        config->file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
        strjoin(config->dataset, "/processed_data", outdat);
        config->lengthmz = mh5getfilelength(config->file_id, outdat);
        inp->dataMZ = calloc(config->lengthmz, sizeof(float));
        inp->dataInt = calloc(config->lengthmz, sizeof(float));
        mh5readfile2d(config->file_id, outdat, config->lengthmz, inp->dataMZ, inp->dataInt);
        printf("Length of Data: %d \n", config->lengthmz);

        //Check the length of the mfile and then read it in.
        if (config->mflag == 1)
        {
            config->mfilelen = mh5getfilelength(config->file_id, "/config/masslist");
            printf("Length of mfile: %d \n", config->mfilelen);
            inp->testmasses = malloc(sizeof(float) * config->mfilelen);
            mh5readfile1d(config->file_id, "/config/masslist", inp->testmasses);
        }
        else {
            inp->testmasses = malloc(sizeof(float) * config->mfilelen);
        }
    }
    else {
        //Calculate the length of the data file automatically
        config->lengthmz = getfilelength(config->infile);
        inp->dataMZ = calloc(config->lengthmz, sizeof(float));
        inp->dataInt = calloc(config->lengthmz, sizeof(float));

        readfile(config->infile, config->lengthmz, inp->dataMZ, inp->dataInt);//load up the data array
        printf("Length of Data: %d \n", config->lengthmz);

        //Check the length of the mfile and then read it in.

        if (config->mflag == 1)
        {
            config->mfilelen = getfilelength(config->mfile);
            printf("Length of mfile: %d \n", config->mfilelen);
        }
        inp->testmasses = malloc(sizeof(float) * config->mfilelen);
        if (config->mflag == 1)
        {
            readmfile(config->mfile, config->mfilelen, inp->testmasses);//read in mass tab
        }
    }


    //This for loop creates a list of charge values
    inp->nztab = calloc(config->numz, sizeof(int));
    for (int i = 0; i < config->numz; i++) { inp->nztab[i] = i + config->startz; }
    //printf("nzstart %d\n",inp.nztab[0]);

    //Test to make sure no charge state is zero
    for (int j = 0; j < config->numz; j++)
    {
        if (inp->nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
    }
    //Test to make sure no two data points has the same x value
    for (int i = 0; i < config->lengthmz - 1; i++)
    {
        if (inp->dataMZ[i] == inp->dataMZ[i + 1]) { printf("Error: Two data points are identical: %f %f \n\n", inp->dataMZ[i], inp->dataMZ[i+1]); exit(104); }
    }

}



void WritePeaks(const Config config, const Decon* decon) {
    char outdat[1024];
    strjoin(config.dataset, "/peaks", outdat);
    float* ptemp = NULL;
    ptemp = calloc(decon->plen * 3, sizeof(float));

    for (int i = 0; i < decon->plen; i++) {
        ptemp[i * 3] = decon->peakx[i];
        ptemp[i * 3 + 1] = decon->peaky[i];
        ptemp[i * 3 + 2] = decon->dscores[i];
    }

    mh5writefile2d_grid(config.file_id, outdat, decon->plen, 3, ptemp);
    free(ptemp);
}


void SetLimits(const Config config, Input* inp)
{
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


Config PostImport(Config config)
{
    //Convert gaussian FWHM to sigma
    if (config.psfun == 0) { config.mzsig = config.mzsig / 2.35482; }
    config.dtsig = config.dtsig / 2.35482;
    //Check whether to turn off or on the config.limitflag. Limit flag on means it will select only the mass values closest to the mass values from mfile.
    if (config.mflag == 1 && config.mtabsig == 0) { config.limitflag = 1; }

    config.numz = config.endz - config.startz + 1;

    //If linflag is active, overwrite speedy}
    if (config.linflag != -1) {
        if (config.linflag != 2) { config.speedyflag = 1; }
        else { config.speedyflag = 0; }
    }
    if (config.speedyflag == 1) { printf("Speedy mode: Assuming Linearized Data\n"); }

    if (config.filetype == 0)
    {
        //Print inputs for check
        printf("infile = %s\n", config.infile);
        //printf("outfile = %s\n", config.outfile);
        //printf("\n");
    }

    //Check to see if the mass axis should be fixed
    if (config.massub < 0 || config.masslb < 0) {
        config.fixedmassaxis = 1;
        config.massub = fabs(config.massub);
        config.masslb = fabs(config.masslb);
    }

    if (config.twaveflag == -1 && config.imflag == 1) { printf("\n\nNeed to define twaveflag for CCS calculation\n\n"); }

    //bug proofing so we don't get a 1/0 problem
    if (config.msig == 0) { config.msig = 0.00001; }
    if (config.zsig == 0) { config.zsig = 0.00001; }
    if (config.massbins == 0) { config.massbins = 1; }

    //Convert aggressiveflag to baselineflag
    if (config.aggressiveflag == 1 || config.aggressiveflag == 2) { config.baselineflag = 1; }
    else { config.baselineflag = 0; }

    //Experimental correction. Not sure why this is necessary.
    if (config.psig < 0) { config.mzsig /= 3; }

    return config;
}



// Unused //

// Manual Assign //


void ManualAssign(float* dataMZ, char* barr, int* nztab, Config config)
{
    unsigned int i, j, k;
    int manlen = 0;
    int lengthmz = config.lengthmz;
    int numz = config.numz;
    if (config.filetype == 1)
    {
        manlen = mh5getfilelength(config.file_id, "/config/manuallist");
    }
    else
    {
        manlen = getfilelength(config.manualfile);
    }
    printf("Length of Manual List: %d \n", manlen);
    float* manualmz = malloc(sizeof(float) * manlen);
    float* manualwin = malloc(sizeof(float) * manlen);
    float* manualassign = malloc(sizeof(float) * manlen);
    if (config.filetype == 1)
    {
        mh5readfile3d(config.file_id, "/config/manuallist", manlen, manualmz, manualwin, manualassign);
    }
    else
    {
        readfile3(config.manualfile, manlen, manualmz, manualwin, manualassign);
    }

    for (i = 0; i < manlen; i++) {
        if (manualassign[i] < 0) { manualmz[i] = -manualmz[i] * 1000.0; }
        //Cheating a bit...make the manualmz very negative if deassign, so that they aren't discovered by the assign loop
    }

    float testmz;
    int closest;
#pragma omp parallel for private (i,j,testmz,closest), schedule(auto)
    for (i = 0;i < lengthmz;i++)
    {
        for (j = 0;j < numz;j++)
        {
            //Manual Assign: The nearest value wins in the case of overlap
            testmz = dataMZ[i];
            closest = nearunsorted(manualmz, testmz, manlen);
            if (fabs(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] == manualassign[closest] && manualassign[closest] > 0)
            {
                barr[index2D(numz, i, j)] = 1;
            }
            else if (fabs(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] != manualassign[closest] && manualassign[closest] > 0)
            {
                barr[index2D(numz, i, j)] = 0;
            }
            //Manual Deassign: Anything within the window is killed
            for (k = 0; k < manlen; k++) {

                if (fabs(manualmz[k] - testmz * (-1000.0)) < manualwin[k] * 1000.0 && (float)nztab[j] == fabs(manualassign[k]) && manualassign[k] < 0)
                {
                    barr[index2D(numz, i, j)] = 0;
                }

            }
        }
    }
    free(manualmz);
    free(manualwin);
    free(manualassign);
    printf("Using Manual Assignments for Some Peaks\n");
}


// Manual Assign //
