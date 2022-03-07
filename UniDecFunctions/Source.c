// File "Source.cpp" 


// you have to include 'pch.h', but also in order for the .dll file to work, you need 
// to go to Properties -> C/C++ -> Precompiled Headers and change the Precompiled Headers setting
// to "Not using precompiled headers." It's a real pain, and it took me like four hours 
// to figure out... 
// Furthermore, you need to change the dllmain.cpp file to dllmain.c file. 

#include "Source.h"
#include "pch.h"
#include <stdio.h>
#include <stdlib.h>

void PrintHelloWorld() {
		if (printf("%s \n", "Hello")) {}; 	
}; 

int ReturnInteger() {
		int i; 
		i = 1;
		return i;
}; 

const char* ReturnCharacterString() {
		// length of C string is +1 what you write because the C compiler 
		// automatically adds the terminating null character at the end of a string
		const char* message = "This string was generated in C and called via the .dll";
		return message;
}; 

char* ReturnModifiedCharacterString(const char* message) {
		// if you're in a C++ file, you need to use memory allocation the C++ way. 
		// E.g. (type)malloc(sizeof(x)); 
		char* modifiedString = (char*)malloc(sizeof(message) + 1);
		char* strReturned; 
		// iterate through the message, changing the vowels. Use a switch for practice. 

		strReturned = modifiedString; 
		return strReturned; 
}; 

float TestArrayOfFloats(float* fArray) {
		float result = 0;
		return result;
}
