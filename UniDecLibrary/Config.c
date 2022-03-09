//
//  Config.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Config.h"

Config SetDefaultConfig()
{
    Config config;
    strcpy(config.infile, "default_data.txt");
    strcpy(config.outfile, "default_output");
    strcpy(config.mfile, "default_mfile.txt");
    config.numit = 50;
    config.endz = 100;
    config.startz = 1;
    config.zsig = 1;
    config.psig = 1;
    config.beta = 0;
    config.mzsig = 15;
    config.msig = 0;
    config.molig = 0;
    config.massub = 5000000;
    config.masslb = 100;
    config.psfun = 0;
    config.mtabsig = 0;
    config.mflag = 0;
    config.massbins = 100;
    config.limitflag = 0;
    config.psthresh = 6;
    config.speedyflag = 0;
    config.aggressiveflag = 0;
    config.adductmass = 1.007276467;
    config.rawflag = 1;
    config.nativezub = 100;
    config.nativezlb = -200;
    config.poolflag = 1;
    config.manualflag = 0;
    config.intthresh = 0;
    config.peakshapeinflate = 1;
    config.killmass = 0;
    config.fixedmassaxis = 0;
    config.isotopemode = 0;
    config.filetype = 0;
    config.imflag = 0;
    config.linflag = -1;
    //IM Parameters
    config.dtsig = 0.2;
    config.csig = 1;
    config.ccsub = 20000;
    config.ccslb = -20000;
    config.ccsbins = 100;
    config.temp = 0;
    config.press = 2;
    config.volt = 50;
    config.tcal1 = 0.3293;
    config.tcal2 = 6.3597;
    config.tcal3 = 0;
    config.tcal4 = 0;
    config.twaveflag = -1;
    config.hmass = 4.002602;
    config.to = 0.97;
    config.len = 0.18202;
    config.edc = 1.57;
    config.nativeccsub = 20000;
    config.nativeccslb = -20000;
    config.zout = 0;
    config.baselineflag = 1;
    config.noiseflag = 0;
    config.metamode = -2;
    config.minmz = -1;
    config.maxmz = -1;
    config.mzbins = 0;
    config.bsub = 0;
    config.datareduction = 0;
    config.peakwin = 500;
    config.peakthresh = 0.1;
    config.exchoice = 0;
    config.exchoicez = 1;
    config.exthresh = 10;
    config.exnorm = 0;
    config.exnormz = 0;
    config.peaknorm = 1;
    config.exwindow = 0;
    config.orbimode = 0;
    config.datanorm = 1;
    //Experimental
    config.filterwidth = 20;
    config.zerolog = -12;
    config.lengthmz = 0;
    config.mfilelen = 0;
    config.isolength = 0;
    // DoubleDec
    config.doubledec = 0; // Consider initializing kernel as well?
    return config;
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

Config LoadConfig(Config config, const char* filename)
{
    // We assume argv[1] is a filename to open
    FILE* file = fopen(filename, "r");

    if (file == 0)
    {
        printf("\nCould not open configuration file: \n \nSomething is wrong...\n\n");
    }
    //Read parameters from configuration file
    //Configuration file should be formatted: name value
    else
    {
        char x[500];
        char y[500];
        //printf("\nRead from file:");
        while (fscanf(file, "%s %500[^\n]", x, y) != EOF)
        {
            //printf( "read in: %s %s \n", x,y );
            if (strstr(x, "input") != NULL) { strcpy(config.infile, y); }// printf(" input");
            if (strstr(x, "output") != NULL) { strcpy(config.outfile, y); }//  printf(" output"); }
            if (strstr(x, "mfile") != NULL) { strcpy(config.mfile, y); config.mflag = 1; }// printf(" mfile"); }
            if (strstr(x, "numit") != NULL) { config.numit = atoi(y); }// printf(" numit"); }
            //if (strstr(x, "numz") != NULL){ config.numz = atoi(y); printf(" numz"); }
            if (strstr(x, "startz") != NULL) { config.startz = atoi(y); }// printf(" startz"); }
            if (strstr(x, "endz") != NULL) { config.endz = atoi(y); }// printf(" endz"); }
            if (strstr(x, "zzsig") != NULL) { config.zsig = atof(y); }// printf(" zzsig"); }
            if (strstr(x, "psig") != NULL) { config.psig = atof(y); }// printf(" psig"); }
            if (strstr(x, "beta") != NULL) { config.beta = atof(y); }// printf(" beta"); }
            if (strstr(x, "mzsig") != NULL) { config.mzsig = atof(y); }// printf(" mzsig"); }
            if (strstr(x, "msig") != NULL) { config.msig = atof(y); }// printf(" msig"); }
            if (strstr(x, "molig") != NULL) { config.molig = atof(y); }// printf(" molig"); }
            if (strstr(x, "massub") != NULL) { config.massub = atof(y); }// printf(" massub"); }
            if (strstr(x, "masslb") != NULL) { config.masslb = atof(y); }// printf(" masslb"); }
            if (strstr(x, "psfun") != NULL) { config.psfun = atoi(y); }// printf(" psfun"); }
            if (strstr(x, "mtabsig") != NULL) { config.mtabsig = atof(y); }// printf(" mtabsig"); }
            if (strstr(x, "massbins") != NULL) { config.massbins = atof(y); }// printf(" massbins"); }
            if (strstr(x, "psthresh") != NULL) { config.psthresh = atof(y); }// printf(" psthresh"); }
            if (strstr(x, "speedy") != NULL) { config.speedyflag = atoi(y); }//  printf(" speedy"); }
            if (strstr(x, "aggressive") != NULL) { config.aggressiveflag = atoi(y); }//  printf(" aggressive"); }
            if (strstr(x, "adductmass") != NULL) { config.adductmass = atof(y); }//  printf(" adductmass"); }
            if (strstr(x, "rawflag") != NULL) { config.rawflag = atoi(y); }// printf(" rawflag"); }
            if (strstr(x, "nativezub") != NULL) { config.nativezub = atof(y); }// printf(" nativezub"); }
            if (strstr(x, "nativezlb") != NULL) { config.nativezlb = atof(y); }// printf(" nativezlb"); }
            if (strstr(x, "poolflag") != NULL) { config.poolflag = atoi(y); }// printf(" poolflag"); }
            if (strstr(x, "manualfile") != NULL) { config.manualflag = 1; strcpy(config.manualfile, y); }//  printf(" manualfile"); }
            if (strstr(x, "intthresh") != NULL) { config.intthresh = atof(y); }// printf(" intthresh"); }
            if (strstr(x, "peakshapeinflate") != NULL) { config.peakshapeinflate = atof(y); }// printf(" peakshapeinflate"); }
            if (strstr(x, "killmass") != NULL) { config.killmass = atof(y); }// printf(" killmass"); }
            if (strstr(x, "isotopemode") != NULL) { config.isotopemode = atoi(y); }// printf(" isotopemode"); }
            if (strstr(x, "orbimode") != NULL) { config.orbimode = atoi(y); }// printf(" orbimode"); }
            if (strstr(x, "imflag") != NULL) { config.imflag = atoi(y); }// printf(" imflag"); }
            if (strstr(x, "linflag") != NULL) { config.linflag = atoi(y); }// printf(" linflag"); }
            //IM Parameters
            if (strstr(x, "csig") != NULL) { config.csig = atof(y); }// printf(" csig"); }
            if (strstr(x, "dtsig") != NULL) { config.dtsig = atof(y); }// printf(" dtsig"); }
            if (strstr(x, "ccsub") != NULL) { config.ccsub = atof(y); }// printf(" ccsub"); }
            if (strstr(x, "ccslb") != NULL) { config.ccslb = atof(y); }// printf(" ccslb"); }
            if (strstr(x, "ccsbins") != NULL) { config.ccsbins = atof(y); }// printf(" ccsbins"); }
            if (strstr(x, "temp") != NULL) { config.temp = atof(y); }// printf(" temp"); }
            if (strstr(x, "pressure") != NULL) { config.press = atof(y); }// printf(" pressure"); }
            if (strstr(x, "volt") != NULL) { config.volt = atof(y); }// printf(" volt"); }
            if (strstr(x, "gasmass") != NULL) { config.hmass = atof(y); }// printf(" gasmass"); }
            if (strstr(x, "tnaught") != NULL) { config.to = atof(y); }// printf(" to"); }
            if (strstr(x, "tcal1") != NULL) { config.tcal1 = atof(y); }// printf(" tcal1"); }
            if (strstr(x, "tcal2") != NULL) { config.tcal2 = atof(y); }// printf(" tcal2"); }
            if (strstr(x, "tcal3") != NULL) { config.tcal3 = atof(y); }// printf(" tcal3"); }
            if (strstr(x, "tcal4") != NULL) { config.tcal4 = atof(y); }// printf(" tcal4"); }
            if (strstr(x, "edc") != NULL) { config.edc = atof(y); }// printf(" edc"); }
            if (strstr(x, "zout") != NULL) { config.zout = atoi(y); }// printf(" zout"); }
            if (strstr(x, "twaveflag") != NULL) { config.twaveflag = atoi(y); }// printf(" twaveflag"); }
            if (strstr(x, "ubnativeccs") != NULL) { config.nativeccsub = atof(y); }// printf(" ubnativeccs"); }
            if (strstr(x, "lbnativeccs") != NULL) { config.nativeccslb = atof(y); }// printf(" lbnativeccs"); }
            if (strstr(x, "driftlength") != NULL) { config.len = atof(y); }// printf(" driftlength"); }
            if (strstr(x, "baselineflag") != NULL) { config.baselineflag = atoi(y); }// printf(" baselineflag"); }
            if (strstr(x, "noiseflag") != NULL) { config.noiseflag = atoi(y); }// printf(" noiseflag"); }
            //Experimental
            if (strstr(x, "filterwidth") != NULL) { config.filterwidth = atoi(y); }// printf(" filterwidth"); }
            if (strstr(x, "zerolog") != NULL) { config.zerolog = atof(y); }// printf(" zerolog"); }
            //Peak Parameters
            if (strstr(x, "peakwindow") != NULL) { config.peakwin = atof(y); }// printf(" peakwindow"); }
            if (strstr(x, "peakthresh") != NULL) { config.peakthresh = atof(y); }// printf(" peakthresh"); }
            if (strstr(x, "peaknorm") != NULL) { config.peaknorm = atoi(y); }// printf(" peaknorm"); }
            // DoubleDec Parameters
            if (strstr(x, "doubledec") != NULL) { config.doubledec = atoi(y); }
            if (strstr(x, "kernel") != NULL) { strcpy(config.kernel, y); }
        }
        //printf("\n\n");
    }
    fclose(file);

    config = PostImport(config);

    return config;

}
