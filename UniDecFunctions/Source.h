#pragma once
/*
* Here's the thing: If you're reading this, you are about to undertake a several hour 
* personal journey into how much frustration can you take before you throw your mouse 
* into the monitor and subsequently throw that into the room. So here's some tips to get you started: 
* 
* 1. Don't try to integrate C++ into C#. It's not worth it. Just change your file extensions
* from .cpp to .c before you even do anything, and you'll save yourself MANY hours of frustration. 
* 2. Resource on default marshaling of types from C to C# https://docs.microsoft.com/en-us/dotnet/standard/native-interop/type-marshaling
* 3. Customiazed structure marshalling for passing complex data back and forth: https://docs.microsoft.com/en-us/dotnet/standard/native-interop/customize-struct-marshaling
* 4. Marshalling different types of arrays: https://docs.microsoft.com/en-us/dotnet/framework/interop/marshaling-different-types-of-arrays
* 5. Oh yeah, here's a resource on marshalling in general: https://docs.microsoft.com/en-us/dotnet/standard/native-interop/type-marshaling
* 6. .lib (static library) vs .dll (dynamically linked library) https://stackoverflow.com/questions/140061/when-to-use-dynamic-vs-static-libraries
* 7. Adding a .dll to a C# project: https://stackoverflow.com/questions/7685718/how-to-add-dll-in-c-sharp-project
* 8. Be prepared to spend a long time suffering. Would highly recommend taking this as a "do some homework style approach." 
* 9. In case you don't take my advice and decided to work with C++ files for some reason: https://www.geeksforgeeks.org/extern-c-in-c/#:~:text=Since%20C%2B%2B%20supports%20function%20overloading,the%20extern%20%E2%80%9CC%E2%80%9D%20block.
* 10. Further C++ help: source on using the header macros: 
 https://stackoverflow.com/questions/14980649/macro-for-dllexport-dllimport-switch
 Summary, find the x_EXPORTS definition in Properties C/C++ -> PreProcessor -> 
 Preprocessor Definitions. You should see something of the form, in this case, 
 UNIDECFUNCTIONS_EXPORTS. The macro below checks if UNIDECFUNCTIONS_EXPORTS is defined, 
 and says that if it is, assign DLLEXPORT the code that signifies that a function is supposed 
 to be exported (__declsepc(dllexport)) 
 If UNIDECFUNCTION_EXPORTS is not defined, then it sets DLLEXPORT to __declspec(dllimport), 
 which allows calling code to import the functions defined from a DLL. 

*/


#include <stdio.h>
#include <stdlib.h>

void PrintHelloWorld(); 
const int ReturnInteger(); 
const char* ReturnModifiedCharacterString(const char* message);
float TestArrayOfFloats(float* fArray); 

