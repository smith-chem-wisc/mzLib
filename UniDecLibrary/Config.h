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

__declspec(dllexport) typedef struct Config Config;

struct Config
{
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
};

// Load config needs to be rewritten so that I can just pass the config struct from the C# calling code.
__declspec(dllexport) void PostImport(Config *config); 
__declspec(dllexport) Config ReturnModifiedConfigToCS(Config config); 
__declspec(dllexport) Config ModifyConfigToDefault(Config* config); 
#endif /* Config_h */
