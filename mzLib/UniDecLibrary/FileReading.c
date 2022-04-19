//
//  FileReading.c
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#include "FileReading.h"

// Reading files
//Reads in x y file.
void readfile(char* infile, int lengthmz, float* dataMZ, float* dataInt)
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
            dataMZ[i] = atof(x);
            dataInt[i] = atof(y);
        }
    }
    fclose(file_ptr);

}

//Reads in x y z file.
void readfile3(char* infile, int lengthmz, float* array1, float* array2, float* array3)
{
    FILE* file_ptr;
    int i;
    char x[500];
    char y[500];
    char z[500];

    file_ptr = fopen(infile, "r");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(10);
    }
    else {

        for (i = 0;i < lengthmz;i++)
        {
            fscanf(file_ptr, "%s %s %s", x, y, z);
            array1[i] = atof(x);
            array2[i] = atof(y);
            array3[i] = atof(z);
        }
    }
    fclose(file_ptr);

}

//Reads in x y z file.
void readfile3bin(char* infile, int lengthmz, float* array1, float* array2, float* array3)
{
    FILE* file_ptr;
    int i;
    float* data;
    data = calloc(lengthmz * 3, sizeof(float));

    file_ptr = fopen(infile, "rb");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(10);
    }
    else {
        fread(data, sizeof(float), lengthmz * 3, file_ptr);
        for (i = 0; i < lengthmz; i++)
        {
            array1[i] = data[i*3];
            array2[i] = data[i*3+1];
            array3[i] = data[i*3+2];
            //printf("%f %f %f\n", array1[i], array2[i], array3[i]);
        }
    }
    fclose(file_ptr);
}


//Reads in single list of values
void readmfile(char* infile, int mfilelen, float* testmasses)
{
    FILE* file_ptr;
    int i;
    char input[500];

    file_ptr = fopen(infile, "r");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(20);
    }
    else {

        for (i = 0;i < mfilelen;i++)
        {
            fscanf(file_ptr, "%s", input);
            testmasses[i] = atof(input);
        }
    }
    fclose(file_ptr);
}

// count the number of lines we have in datafile
int getfilelength(const char* infile)
{
    FILE* file_ptr;
    int l = 0;
    char input[501];
    file_ptr = fopen(infile, "r");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(9);
    }
    else
    {
        while (fscanf(file_ptr, "%500[^\n]%*c", input) != EOF)
        {
            l += 1;
        }
    }
    fclose(file_ptr);
    return l;
}

// count the number of lines we have in datafile
int getfilelengthbin(const char* infile, const int size, const int width)
{
    FILE* file_ptr;
    int l = 0;
    char input[501];
    file_ptr = fopen(infile, "rb");

    if (file_ptr == 0)
    {
        printf("Could not open %s file\n", infile);
        exit(9);
    }
    else
    {
        
        fseek(file_ptr, 0, SEEK_END);
        l = ftell(file_ptr);
        l /= (size * width);
    }
    fclose(file_ptr);
    return l;
}


