#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "UniDec.h"
#include "UD_score.h"
#include "Config.h"
#include "Decon.h"
#include "Input.h"
#include "TestMasses.h"
#include "MathUtilities.h"
#include "Convolution.h"
#include "ArgMax.h"
#include "PointSmoothing.h"
#include "ErrorFunctions.h"
#include "Decon.h"
#include "Convolution.h"
#include "MZPeak.h"
#include "Transforms.h"
#include "Interpolation.h"
#include "Isotopes.h"
#include "FitFunctions.h"
#include "Normalization.h"
#include "UD_analysis.h"

__declspec(dllexport) int DeclarePointersAndFree(void); 