/*
* UniDec_Main.h
*
*  Created on : 29 December 2016
* Author : Michael.Marty
*/

//
// Copyright 2015 University of Oxford
// Copyright 2016 University of Arizona
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "UniDec.h"
#include "UD_H5_IO.h"
#include "UD_score.h"


Decon MainDeconvolution(const Config config, const Input inp, const int silent, const int verbose);
int run_unidec(int argc, char *argv[], Config config);
void RunAutotune(Config *config, Input *inp, Decon *decon);
Decon MainDeconvolution(const Config config, const Input inp, const int silent, const int verbose);
