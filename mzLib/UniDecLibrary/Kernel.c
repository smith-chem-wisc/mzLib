//
//  Kernel.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "Kernel.h"

// Kernel //

//Reads in kernel file.
void readkernel(const char* infile, int lengthmz, double* datax, double* datay)
{
    FILE* file_ptr;
    int i;
    char x[500];
    char y[500];

    file_ptr = fopen(infile, "r");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(10);
    }
    else {

        for (i = 0;i < lengthmz;i++)
        {
            fscanf(file_ptr, "%s %s", x, y);
            datax[i] = atof(x);
            datay[i] = atof(y);
        }
    }
    fclose(file_ptr);

}


// Kernel //

