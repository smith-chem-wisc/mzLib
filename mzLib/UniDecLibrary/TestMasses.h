//
//  TestMasses.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef TestMasses_h
#define TestMasses_h

#include <stdio.h>
void KillMass(float killmass, int lengthmz, int numz, char* barr, int* nztab, float adductmass, float* dataMZ, float psfun, float mzsig);

void TestMassListWindowed(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab, float* testmasses, int mfilelen, float mtabsig);
void TestMassListLimit(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab, int* testmasspos, int mfilelen);
void TestMass(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab);

#endif /* TestMasses_h */
