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
#include "UniDec_Main.h"

__declspec(dllexport) int MemoryAllocationOfBarr(void); 
__declspec(dllexport) int AllocateMemoryToPointersThenFree(); 
__declspec(dllexport) char UseMemcpyInC();
__declspec(dllexport) char UseMemcpyWithInpAndConfigObjects(); 
__declspec(dllexport) int MemoryObjectAllocationToHeap(Config config, Input inp); 
__declspec(dllexport) int TestSetStartEnds(Input inp, Config config); 
__declspec(dllexport) int MemoryObjectAllocationToHeapConfigPtr(Config* config, Input inp); 
__declspec(dllexport) int TestFreeDecon(); 
__declspec(dllexport) int TestSetupAndAllocateMemoryToDecon(); 
__declspec(dllexport) Decon TestSetupAndReturnDecon(); 
__declspec(dllexport) int MainDeconWithMinimalControlFlow(Config config, Input inp); 
__declspec(dllexport) Decon RunUniDecWithTestMainDeconAlgo(Input inp, Config config);