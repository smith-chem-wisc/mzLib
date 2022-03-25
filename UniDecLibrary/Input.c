//
//  Input.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include <stdio.h>
#include "Input.h"
#include <string.h>
#include <stdlib.h>


Input SetupInputs(void){
        Input inp;
        inp.dataMZ = NULL;
        inp.dataInt = NULL;
        inp.testmasses = NULL;
        inp.nztab = NULL;
        inp.mtab = NULL;
        inp.barr = NULL;
        inp.isotopepos = NULL;
        inp.isotopeval = NULL;
        //Set default isotope parameters. These can be overwritten by config file.
    
        return inp;
}

void FreeInputs(Input inp){
    free(inp.dataMZ);
    free(inp.dataInt);
    free(inp.nztab);
    free(inp.mtab);
    free(inp.testmasses);
    free(inp.barr);
    free(inp.isotopepos);
    free(inp.isotopeval);
}
Input SetupInputsSafe(Input* inp) {
    inp->dataMZ = NULL;
    inp->dataInt = NULL;
    inp->testmasses = NULL;
    inp->nztab = NULL;
    inp->mtab = NULL;
    inp->barr = NULL;
    inp->isotopepos = NULL;
    inp->isotopeval = NULL;
    //Set default isotope parameters. These can be overwritten by config file.
    float temp[10] = { 1.00840852e+00, 1.25318718e-03, 2.37226341e+00, 8.19178000e-04, -4.37741951e-01, 6.64992972e-04, 9.94230511e-01, 4.64975237e-01, 1.00529041e-02, 5.81240792e-01 };
    memcpy(inp->isoparams, temp, sizeof(temp));
    return *inp;
} 


