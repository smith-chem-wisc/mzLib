//
//  Config.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Config_h
#define Config_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct Config Config;

struct Config
{
    char infile[500]; // remove infile
    char outfile[500]; // remove outfile
    int numit;// number iterations
    int numz; // number of charge states
    int endz; // final charge state
    int startz; // initial charge state
    float zsig; //
    float psig;
    float beta;
    float mzsig;
    float msig;
    float molig;
    float massub;
    float masslb;
    int psfun;
    float mtabsig;
    char mfile[500]; // remove mfile
    char manualfile[500]; // remove manualfile
    int mflag;
    float massbins;
    int limitflag;
    float psthresh;
    int speedyflag;
    int linflag;
    int aggressiveflag;
    float adductmass;
    int rawflag;
    float nativezub;
    float nativezlb;
    int poolflag;
    int manualflag;
    float intthresh;
    float peakshapeinflate;
    float killmass;
    int fixedmassaxis;
    int isotopemode;
    int filetype;
    int imflag;
    //IM Parameters
    float dtsig;
    float csig;
    float ccsub;
    float ccslb;
    float ccsbins;
    float temp;
    float press;
    float volt;
    float tcal1;
    float tcal2;
    float tcal3;
    float tcal4;
    float twaveflag;
    float hmass;
    float to;
    float len;
    float edc;
    float nativeccsub;
    float nativeccslb;
    int baselineflag;
    int noiseflag;
    int zout;
    int metamode;
    float minmz;
    float maxmz;
    int mzbins;
    float bsub;
    float datareduction;
    float peakwin;
    float peakthresh;
    float exwindow;
    int exchoice;
    int exchoicez;
    float exthresh;
    int exnorm;
    int exnormz;
    int peaknorm;
    int orbimode;
    int datanorm;
    //Experimental Parameters
    int filterwidth;
    float zerolog;
    int lengthmz;
    int mfilelen;
    int isolength;
    // DoubleDec Parameters
    int doubledec;
    char kernel[500];
    char dataset[1024];
};

Config SetDefaultConfig(void);
// Load config needs to be rewritten so that I can just pass the config struct from the C# calling code.
Config LoadConfig(Config config, const char* filename);

#endif /* Config_h */
