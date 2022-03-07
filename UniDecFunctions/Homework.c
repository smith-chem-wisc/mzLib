#include "Homework.h"
#include "pch.h"
// If you create a pointer in a function, you also need to write the code 
// to delete that pointer. 
int AddTwoIntegers(int i, int j) {
	return i + j; 
}

double MultiplyTwoDouble(double i, double j) {
	return i * j;
}

char* EchoMessage(char* message) {
	return message; 
}

// in C compilation int* is equivalent to int[]. 
int* PassIntArrayModifyThenReturnArray(int* ArrayOfInts, int length) {
	for (int i = 0; i < length; i++) {
		ArrayOfInts[i] *= 2; 
	}
	// note that we are never reallocating new memory, all we're doing is 
	// modifying the memory specified by the pointer. If we need to create a new memory location,
	// we will need to also write a method that de-allocates that memory location and call that from the 
	// C# code as well. This would be annoying to do, and very easy to forget, so please try to stick to only 
	// initializing object in the C#.  
	return ArrayOfInts; 
}
TestStruct InitTestStruct() {
	int val1 = 1; 

	TestStruct test; 
	test.val1 = val1; 
	return test; 
}

TestStruct InitializeAndReturnStruct() {
	return InitTestStruct(); 
}

char* SendingAMessage() {
	char message[100] = "It's not about the money... It's about sending a message."; 
	return &message; 
}

TestStruct2 InitTestStruct2(int val1, int* array1) {
	TestStruct2 testStruct2; 
	testStruct2.val1 = val1; 
	testStruct2.array1 = array1; 
	return testStruct2; 
}
