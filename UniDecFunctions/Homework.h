#pragma once
/*
* The goal of this project is to get you to work out how to pass data backand forth between
* unmanaged C code and managed C# code. 
* 
* You really need to think about type marshalling in the C# code. 
*/ 

// To get this function to run: Build the C .dll in Release mode. You need to add __declspec(dllexport) to 
// all function you are planning to export. 

#include <stdio.h>
#include <stdlib.h>

__declspec(dllexport) int AddTwoIntegers(int i, int j);
__declspec(dllexport) double MultiplyTwoDouble(double i, double j);
__declspec(dllexport) char* EchoMessage(char* message);
__declspec(dllexport) const char* EchoMessageConst(const char* message);
// Pass a pointer to a float[] from C# to C, then multiply 
// each element in the array by 2. Then return the pointer to the 
// now modified array, and fill a new float[] with the modified values. 
// You'll need to figure out how to pass by reference in C#: 
// https://stackoverflow.com/questions/58564167/c-sharp-pinvoke-pass-reference-to-int
__declspec(dllexport) int* PassIntArrayModifyThenReturnArray(int* ArrayOfInts, int length); 

__declspec(dllexport) typedef struct TestStruct {
	int val1; 
} TestStruct;

TestStruct InitTestStruct();

__declspec(dllexport) TestStruct InitializeAndReturnStruct(); 

__declspec(dllexport) char* SendingAMessage(); 

__declspec(dllexport) typedef struct TestStruct2 {
	int val1;
	int* array1; 
} TestStruct2;

__declspec(dllexport) TestStruct2 InitTestStruct2(int val1, int* array1); 

