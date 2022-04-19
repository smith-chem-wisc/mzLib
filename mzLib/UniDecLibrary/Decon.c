//
//  Decon.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Decon.h"

Decon SetupDecon(void){
    Decon decon;
    decon.fitdat = NULL;
    decon.baseline = NULL;
    decon.noise = NULL;
    decon.massgrid = NULL;
    decon.massaxis = NULL;
    decon.massaxisval = NULL;
    decon.blur = NULL;
    decon.newblur = NULL;
    decon.peakx = NULL;
    decon.peaky = NULL;
    decon.dscores = NULL;
    decon.error = 0;
    decon.rsquared = 0;
    decon.iterations = 0;
    decon.uniscore = 0;
    decon.conv = 0;
    decon.threshold = 0;
    decon.mlen = 0;
    decon.plen = 0;
    return decon;
}

void FreeDecon(Decon decon){
    free(decon.fitdat);
    free(decon.baseline);
    free(decon.noise);
    free(decon.massgrid);
    free(decon.massaxis);
    free(decon.massaxisval);
    free(decon.blur);
    free(decon.newblur);
    free(decon.peakx);
    free(decon.peaky);
    free(decon.dscores);
}


