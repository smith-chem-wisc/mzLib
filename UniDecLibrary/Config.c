//
//  Config.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Config.h"

Config SetDefaultConfig(Config config)
{
    config.numit = 50;
    config.numz = 1; 
    config.endz = 100;
    config.startz = 1;
    config.zsig = (float)1;
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
Config ModifyConfigToDefault(Config* config) {
    config->numit = 50;
    config->numz = 1;
    config->endz = 100;
    config->startz = 1;
    config->zsig = (float)1;
    config->psig = 1;
    config->beta = 0;
    config->mzsig = 15;
    config->msig = 0;
    config->molig = 0;
    config->massub = 5000000;
    config->masslb = 100;
    config->psfun = 0;
    config->mtabsig = 0;
    config->mflag = 0;
    config->massbins = 100;
    config->limitflag = 0;
    config->psthresh = 6;
    config->speedyflag = 0;
    config->aggressiveflag = 0;
    config->adductmass = 1.007276467;
    config->rawflag = 1;
    config->nativezub = 100;
    config->nativezlb = -200;
    config->poolflag = 1;
    config->manualflag = 0;
    config->intthresh = 0;
    config->peakshapeinflate = 1;
    config->killmass = 0;
    config->fixedmassaxis = 0;
    config->isotopemode = 0;
    config->imflag = 0;
    config->linflag = -1;
    //IM Prameters
    config->dtsig = 0.2;
    config->csig = 1;
    config->ccsub = 20000;
    config->ccslb = -20000;
    config->ccsbins = 100;
    config->temp = 0;
    config->press = 2;
    config->volt = 50;
    config->tcal1 = 0.3293;
    config->tcal2 = 6.3597;
    config->tcal3 = 0;
    config->tcal4 = 0;
    config->twaveflag = -1;
    config->hmass = 4.002602;
    config->to = 0.97;
    config->len = 0.18202;
    config->edc = 1.57;
    config->nativeccsub = 20000;
    config->nativeccslb = -20000;
    config->zout = 0;
    config->baselineflag = 1;
    config->noiseflag = 0;
    config->metamode = -2;
    config->minmz = -1;
    config->maxmz = -1;
    config->mzbins = 0;
    config->bsub = 0;
    config->datareduction = 0;
    config->peakwin = 500;
    config->peakthresh = 0.1;
    config->exchoice = 0;
    config->exchoicez = 1;
    config->exthresh = 10;
    config->exnorm = 0;
    config->exnormz = 0;
    config->peaknorm = 1;
    config->exwindow = 0;
    config->orbimode = 0;
    config->datanorm = 1;
    //Expeimental
    config->filterwidth = 20;
    config->zerolog = -12;
    config->lengthmz = 0;
    config->mfilelen = 0;
    config->isolength = 0;
    // DouleDec
    config->doubledec = 0; // Consider initializing kernel as well?    
    return *config; 
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

Config ReturnModifiedConfigToCS(Config config) {
    Config* ptr = &config; 
    ptr->adductmass = (float)23.0; 
    return config; 
}
