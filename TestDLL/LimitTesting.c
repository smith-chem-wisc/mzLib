void TestMassListWindowed(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab, float* testmasses, int mfilelen, float mtabsig)
{
	unsigned int i, j;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			barr[index2D(numz, i, j)] = 0;
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub && testmass>masslb && nztab[j]<nativelimit + nativezub && nztab[j]>nativelimit + nativezlb)
			{
				if (neartest(testmasses, testmass, mfilelen, mtabsig) == 1)
				{
					barr[index2D(numz, i, j)] = 1;
				}
			}
		}
	}
}

void TestMassListLimit(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab, int* testmasspos, int mfilelen)
{
	unsigned int i, j, k;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			barr[index2D(numz, i, j)] = 0;
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub && testmass>masslb && nztab[j]<nativelimit + nativezub && nztab[j]>nativelimit + nativezlb)
			{
				for (k = 0; k < mfilelen; k++)
				{
					if (testmasspos[index2D(numz, k, j)] == i)
					{
						barr[index2D(numz, i, j)] = 1;
					}
				}
			}
		}
	}
}

void TestMass(int lengthmz, int numz, char* barr, float* mtab, float nativezub, float nativezlb, float massub, float masslb, int* nztab)
{
	unsigned int i, j;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub && testmass>masslb && nztab[j]<nativelimit + nativezub && nztab[j]>nativelimit + nativezlb)
			{
				barr[index2D(numz, i, j)] = 1;
			}
			else
			{
				barr[index2D(numz, i, j)] = 0;
			}
		}
	}
}